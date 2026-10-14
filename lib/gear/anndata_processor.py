"""
AnndataProcessor - Process expression datasets in various formats to H5AD.

Handles conversion, sanitization, and primary analysis of expression datasets
in formats: H5AD, 3-tab, Excel, MEX.
"""

import gc
import json
import os
import sys
import tarfile
import zipfile
from pathlib import Path

import anndata
import geardb
import pandas as pd
import scanpy as sc
from gear.primary_analysis import (
    PrimaryAnalysisProcessingError,
    add_primary_analysis_to_dataset,
)
from gear.utils.gene_mapping import (
    map_gene_symbols_via_mygene,
    update_var_with_ensembl_ids,
)
from gear.utils.job_coordination import UploadCancelledError, raise_if_upload_deleted
from gear.utils.archives import ARCHIVE_READ_ERRORS, archive_error_message
from gear.utils.obs import (
    flag_ambiguous_obs_columns,
    standardize_and_sanitize_obs,
)
from scipy import sparse

_this_dir = Path(__file__).resolve().parent
_root_dir = _this_dir.parents[1]
UPLOADS_BASE_DIR = _root_dir / "www" / "uploads" / "files"


def write_status(status_file, status):
    """Write status dictionary to JSON file."""
    with open(status_file, 'w') as f:
        f.write(json.dumps(status, indent=2))

def process_anndata_synchronously(job_id, share_uid, staging_area, status_file, dataset_uid, dataset_format, perform_primary_analysis):
    """Process anndata upload synchronously."""
    processor = AnndataProcessor(
        job_id=job_id,
        share_uid=share_uid,
        staging_area=staging_area,
        status_file=status_file,
        dataset_uid=dataset_uid,
    )
    return processor.process(dataset_format, perform_primary_analysis)

def clean_chunk(chunk: pd.DataFrame) -> pd.DataFrame:
    """
    Clean a dataframe chunk by removing whitespace and converting to numeric.

    Args:
        chunk: DataFrame chunk to clean

    Returns:
        Cleaned DataFrame
    """
    chunk = chunk.replace(r'^\s+|\s+$', '', regex=True)
    chunk = chunk.apply(
        lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x)
    )
    chunk = chunk.apply(pd.to_numeric, errors='coerce').fillna(0)
    return chunk


def package_content_type(filenames: list[str]) -> str | None:
        #print("DEBUG: filenames", file=sys.stderr, flush=True)
        #print(filenames, file=sys.stderr, flush=True)
        """
        mex:
        matrix.mtx
        barcodes.tsv
        genes.tsv
        (or matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)

        threetab:
        expression.tab
        genes.tab
        observations.tab

        None is returned if neither of these is true

        Added NEMO file format functionality.
        DataMTX.tab -> expression.tab
        COLmeta.tab -> observations.tab
        ROWmeta.tab -> genes.tab
        """
        if 'expression.tab' in filenames and 'genes.tab' in filenames and 'observations.tab' in filenames:
            return 'threetab'

        # MEX files may sit inside a folder in the archive (e.g. filtered_feature_bc_matrix/),
        #  and may be the legacy uncompressed set or the gzipped Cell Ranger v3+ set
        basenames = {Path(f).name for f in filenames}
        if {'matrix.mtx', 'barcodes.tsv', 'genes.tsv'} <= basenames:
            return 'mex'
        if {'matrix.mtx.gz', 'barcodes.tsv.gz', 'features.tsv.gz'} <= basenames:
            return 'mex'

        if 'DataMTX.tab' in filenames and 'COLmeta.tab' in filenames and 'ROWmeta.tab' in filenames:
            return 'threetab'

        return None

class ProcessingError(Exception):
    """Raised when dataset processing fails."""
    pass

class AnndataProcessor:
    """Process and convert expression datasets to H5AD format."""

    def __init__(
        self,
        job_id: str,
        share_uid: str,
        staging_area: Path,
        status_file: Path,
        dataset_uid: str,
    ) -> None:
        """
        Initialize the processor.

        Args:
            job_id: Unique job identifier
            share_uid: Share UID for the dataset
            staging_area: Directory containing uploaded files
            status_file: Path to status.json for progress updates
            dataset_uid: Dataset UID for primary analysis
        """
        resolved_staging_area = Path(staging_area).resolve()
        if not resolved_staging_area.is_relative_to(UPLOADS_BASE_DIR.resolve()):
            raise ProcessingError(f"Invalid staging area path: {staging_area}")

        self.job_id = job_id
        self.share_uid = share_uid
        self.staging_area = resolved_staging_area
        self.status_file = status_file
        self.dataset_uid = dataset_uid
        self.status = {
            "job_id": self.job_id,
            "status": "processing",
            "message": "",
            "progress": 0,
        }

    def process(
        self,
        dataset_format: str,
        perform_primary_analysis: bool = False,
    ) -> dict:
        """
        Process the uploaded dataset.

        Args:
            dataset_format: Format of the dataset (h5ad, threetab, excel, mex)
            perform_primary_analysis: Whether to perform primary analysis after processing

        Returns:
            Result dictionary with 'success' and 'message' keys
        """
        try:
            self._process_by_format(dataset_format)

            if perform_primary_analysis:
                self._update_progress(66, "Performing primary analysis...")
                try:
                    add_primary_analysis_to_dataset(
                        self.dataset_uid,
                        self.share_uid,
                        self.staging_area,
                        dataset_format,
                    )
                except PrimaryAnalysisProcessingError as e:
                    raise ProcessingError(
                        f"Primary analysis (clustering and dimensionality reduction) failed: {e}. "
                        "This can happen if the dataset has very few cells/genes or other data quality "
                        "issues. If the dataset looks correct to you, please contact the gEAR team for "
                        f"help and reference share ID {self.share_uid}."
                    )

            message =  "Dataset processed successfully."

            self._update_progress(100, message)
            self._update_status("complete", message)
            return {"success": 1, "message": message}

        except UploadCancelledError as e:
            return {"success": 0, "cancelled": True, "message": str(e)}
        except ProcessingError as e:
            # A write into a deleted staging directory can surface as a ProcessingError too
            if not self.staging_area.is_dir():
                return self._cancelled_result()
            self._update_status("error", str(e))
            return {"success": 0, "message": str(e)}
        except Exception as e:
            if not self.staging_area.is_dir():
                return self._cancelled_result()
            import traceback
            traceback.print_exc()
            message = (
                "An unexpected internal error occurred while processing this dataset. "
                "Please contact the gEAR team to resolve this issue (share ID: "
                f"{self.share_uid}). (Technical details: {e})"
            )
            self._update_status("error", message)
            return {"success": 0, "message": message}
        finally:
            gc.collect()

    def _process_by_format(self, dataset_format: str) -> Path:
        """Route to appropriate processor based on dataset format."""
        if dataset_format == "h5ad":
            return self._process_h5ad()
        elif dataset_format == "mex_3tab":
            return self._process_mex_3tab()
        elif dataset_format == "rds":
            return self._process_seurat()
        elif dataset_format == "excel":
            return self._process_excel()
        else:
            # The uploader should have already rejected unsupported formats before
            # queuing this job, so reaching here points to an internal routing bug.
            raise ProcessingError(
                f"Internal error: received an unrecognized dataset format ('{dataset_format}') for "
                f"processing. Please contact the gEAR team to resolve this issue (share ID: "
                f"{self.share_uid})."
            )

    def _process_h5ad(self) -> Path:
        """Process .h5ad file with backed mode for memory efficiency."""
        self._update_progress(5, "Reading H5AD file...")

        filepath = self.staging_area / f"{self.share_uid}.h5ad"

        # Read in backed mode to avoid loading full X matrix into memory
        try:
            adata = anndata.read_h5ad(filepath, backed='r')
        except Exception as e:
            raise ProcessingError(
                f"Could not read the uploaded H5AD file: {e}. Please verify the file is a valid, "
                "uncorrupted AnnData/H5AD file (for example, by confirming it opens with "
                "scanpy.read_h5ad() or anndata.read_h5ad() locally) and re-upload it."
            )

        self._update_progress(15, "Standardizing and flagging observation metadata...")

        # obs/var are metadata only (small); safe to mutate in place on the backed object
        self._sanitize_and_flag_obs_columns(adata)

        if "gene_symbol" not in adata.var.columns:
            self._update_progress(25, "Mapping gene symbols via Ensembl...")
            adata.var = self._update_var_with_ensembl_ids(adata.var)

        self._update_progress(50, "Writing H5AD file...")

        # Write directly from the backed object so anndata streams X from the
        # source file in chunks rather than loading the full matrix into memory
        h5ad_temp = self.staging_area / f"{self.share_uid}.new.h5ad"
        try:
            adata.write(h5ad_temp, compression='gzip')
        except Exception as e:
            raise ProcessingError(
                f"Failed to write the processed H5AD file: {e}. This is an internal issue, not "
                f"something wrong with your data — please contact the gEAR team for help and "
                f"reference share ID {self.share_uid}."
            )
        finally:
            adata.file.close()

        filepath.unlink()
        h5ad_temp.rename(filepath)

        self._update_progress(65, "H5AD processing complete.")
        return filepath

    def _process_mex_3tab(self) -> Path:
        """Extract an uploaded MEX or 3-tab archive (.tar, .tar.gz or .zip) and process its contents."""
        # store_expression_dataset.cgi saves the upload as <share_uid>.<extension>
        filename = None
        for extension in ("tar.gz", "tar", "zip"):
            candidate = self.staging_area / f"{self.share_uid}.{extension}"
            if candidate.exists():
                filename = candidate
                break

        if filename is None:
            raise ProcessingError(
                "The uploaded archive for this dataset could not be found on the server. "
                "This is usually a transient upload issue — please try re-uploading the dataset. "
                f"If it keeps happening, contact the gEAR team and reference share ID {self.share_uid}."
            )

        compression_format = 'zip' if filename.suffix == '.zip' else 'tarball'

        files_extracted = []

        if compression_format == 'tarball':
            try:
                # "r:*" detects compression from the file contents, so both plain and gzipped
                #  tarballs open regardless of the extension the user gave them
                with tarfile.open(filename, "r:*") as tf:
                    for entry in tf:
                        # Extraction recreates missing directories, so stop if the upload was deleted
                        raise_if_upload_deleted(self.staging_area)
                        tf.extract(entry, path=self.staging_area)

                        # Nemo suffixes
                        nemo_suffixes = ['DataMTX.tab', 'COLmeta.tab', 'ROWmeta.tab']
                        suffix_found = None

                        for suffix in nemo_suffixes:
                            if entry.name.endswith(suffix):
                                suffix_found = suffix
                                # Rename the file to the appropriate name
                                old_name = self.staging_area / entry.name
                                new_name = self.staging_area / suffix
                                old_name.replace(new_name)

                        if suffix_found is not None:
                            files_extracted.append(suffix_found)
                        else:
                            files_extracted.append(entry.name)
            except ARCHIVE_READ_ERRORS as e:
                # Includes a truncated .tar.gz, which fails with EOFError rather than ReadError
                raise ProcessingError(archive_error_message(e, str(filename)))

        if compression_format == 'zip':
            try:
                with zipfile.ZipFile(filename) as zf:
                    for entry in zf.infolist():
                        # Extraction recreates missing directories, so stop if the upload was deleted
                        raise_if_upload_deleted(self.staging_area)
                        zf.extract(entry, path=self.staging_area)

                        # Nemo suffixes
                        nemo_suffixes = ['DataMTX.tab', 'COLmeta.tab', 'ROWmeta.tab']
                        suffix_found = None

                        for suffix in nemo_suffixes:
                            if entry.filename.endswith(suffix):
                                suffix_found = suffix
                                # Rename the file to the appropriate name
                                old_name = self.staging_area / entry.filename
                                new_name = self.staging_area / suffix
                                old_name.replace(new_name)

                        if suffix_found is not None:
                            files_extracted.append(suffix_found)
                        else:
                            files_extracted.append(entry.filename)
            except ARCHIVE_READ_ERRORS as e:
                raise ProcessingError(archive_error_message(e, str(filename)))

        # Determine the dataset type
        dataset_type = package_content_type(files_extracted)

        if dataset_type is None:
            raise ProcessingError(
                "Could not determine the dataset format from the files in your archive. For a "
                "3-tab dataset, it must contain expression.tab, genes.tab, and observations.tab "
                "(or the NeMO-style DataMTX.tab/ROWmeta.tab/COLmeta.tab equivalents). For a MEX "
                "dataset, it must contain matrix.mtx, barcodes.tsv, and genes.tsv. Please check "
                "your archive contains one of these complete file sets and re-upload it."
            )

        # Call the appropriate function
        if dataset_type == 'threetab':
            return self._process_threetab()
        elif dataset_type == 'mex':
            return self._process_mex()

        # Unreachable: package_content_type() only ever returns 'threetab', 'mex', or None
        # (handled above), so getting here indicates an internal logic error.
        raise ProcessingError(
            f"Internal error: detected dataset type '{dataset_type}' has no matching processor. "
            f"Please contact the gEAR team for help and reference share ID {self.share_uid}."
        )


    def _process_threetab(self) -> Path:
        """Process 3-tab format (expression, genes, observations) with chunking."""
        self._update_progress(5, "Extracting and validating files...")

        # Extract/find the three required files
        expression_matrix_path, obs, var = self._extract_threetab_files()

        self._update_progress(25, "Processing expression matrix in chunks...")

        # Process expression matrix in chunks
        chunk_size = 500
        expression_matrix = self._read_expression_matrix_chunks(
            expression_matrix_path, chunk_size, total_rows=len(var)
        )

        # Create AnnData object
        adata = anndata.AnnData(X=expression_matrix, obs=var, var=obs)
        adata = adata.transpose()

        self._update_progress(40, "Standardizing and flagging observation metadata...")
        self._sanitize_and_flag_obs_columns(adata)

        self._update_progress(50, "Writing H5AD file...")

        h5ad_path = self.staging_area / f"{self.share_uid}.h5ad"
        adata.write(h5ad_path, compression='gzip')

        self._update_progress(65, "3-tab processing complete.")
        return h5ad_path

    def _extract_threetab_files(self) -> tuple[Path, pd.DataFrame, pd.DataFrame]:
        """Extract and read the three required 3-tab files."""
        expression_matrix_path = None
        obs = None
        var = None

        for infile in self.staging_area.iterdir():
            if infile.name.startswith('.'):
                continue

            filepath = infile

            if infile.name in ['expression.tab', 'DataMTX.tab']:
                expression_matrix_path = filepath
            elif infile.name in ['observations.tab', 'COLmeta.tab']:
                obs = pd.read_table(filepath, sep='\t', index_col=0, header=0)
            elif infile.name in ['genes.tab', 'ROWmeta.tab']:
                var = pd.read_table(filepath, sep='\t', index_col=0, header=0)

        if obs is None:
            raise ProcessingError(
                "No observations file found in your archive (expected observations.tab or "
                "COLmeta.tab). Please include this file and re-upload."
            )
        if var is None:
            raise ProcessingError(
                "No genes file found in your archive (expected genes.tab or ROWmeta.tab). "
                "Please include this file and re-upload."
            )
        if expression_matrix_path is None:
            raise ProcessingError(
                "No expression matrix file found in your archive (expected expression.tab or "
                "DataMTX.tab). Please include this file and re-upload."
            )

        return expression_matrix_path, obs, var

    def _process_seurat(self) -> Path:
        import gear.seuratuploader as SeuratUploader

        if self.share_uid is None:
            raise ProcessingError(
                f"Internal error: no share ID was provided for this upload job. Please contact "
                f"the gEAR team to resolve this issue (job ID: {self.job_id})."
            )

        # Take in an RDS file, convert to anndata, update the obs metadata based on reductions,
        # convert gene symbols to ensemble IDs, and write to an updated h5ad file.
        self._update_progress(5, "Starting Seurat processing...")
        seurat_filepath = self.staging_area / f"{self.share_uid}.rds"

        # seurat to anndata uses rpy2 to convert the RDS to anndata
        # filepath name has "tmp_" appended in front
        try:
            adata_filepath = SeuratUploader.seurat_to_anndata(
                str(seurat_filepath),
                self.share_uid,
                str(self.staging_area),
                progress_callback=self._update_progress,
            )
        except Exception as e:
            raise ProcessingError(
                f"Could not convert the uploaded RDS file to H5AD: {e}. Please verify this is a "
                "valid Seurat object saved with saveRDS() (from a supported Seurat version), or "
                f"contact the gEAR team for help and reference share ID {self.share_uid} if you believe "
                "the file is valid."
            )

        self._update_progress(15, "Reading converted h5ad file...")
        try:
            adata = anndata.read_h5ad(adata_filepath)
        except Exception as e:
            raise ProcessingError(
                f"Internal error: could not read the H5AD file produced from your RDS conversion: "
                f"{e}. Please contact the gEAR team to resolve this issue (share ID: "
                f"{self.share_uid})."
            )

        # Record the resulting shape so a wildly-off conversion (e.g. from a truncated/corrupted
        # upload) is visible in status.json/logs immediately, rather than only discoverable by
        # manually inspecting the output file after the job reports "complete".
        print(
            f"INFO: Converted RDS for share ID {self.share_uid} to AnnData with shape "
            f"{adata.shape[0]} obs x {adata.shape[1]} vars.",
            file=sys.stderr,
        )
        self._update_progress(
            20, f"Converted to {adata.shape[0]} observations x {adata.shape[1]} variables..."
        )

        # Update obs metadata based on reductions
        self._update_progress(25, "Updating metadata from reductions...")
        try:
            adata = SeuratUploader.reduction_to_metadata(adata)
        except Exception as e:
            raise ProcessingError(
                f"Could not read the dimensionality-reduction embeddings (e.g. PCA/UMAP) from "
                f"your Seurat object: {e}. Please verify the object has valid reductions stored, "
                f"or contact the gEAR team for help and reference share ID {self.share_uid}."
            )

        self._update_progress(30, "Standardizing and flagging observation metadata...")
        try:
            self._sanitize_and_flag_obs_columns(adata)
        except Exception as e:
            raise ProcessingError(
                f"Could not sanitize the observation metadata from your Seurat object: {e}. "
                f"Please verify the object has valid obs data, or contact the gEAR team for help "
                f"and reference share ID {self.share_uid}."
            )

        # Convert gene symbols to ensemble IDs
        metadata_file = self.staging_area / 'metadata.json'
        if not metadata_file.is_file():
            raise ProcessingError(
                f"Internal error: no metadata file was found for this upload job. Please contact "
                f"the gEAR team for help and reference share ID {self.share_uid}."
            )
        # get organism_id by converting sample_taxid(needed for some but not all spatial handlers)
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        self._update_progress(35, "Converting gene symbols to Ensembl IDs...")
        sample_taxid = metadata.get("sample_taxid", None)
        try:
            adata = SeuratUploader.genes_to_ensembl(adata, sample_taxid)
            if adata is None:
                raise ProcessingError("genes_to_ensembl Anndata conversion returned None")
        except Exception as e:
            raise ProcessingError(
                f"Could not map gene symbols to Ensembl IDs: {e}. Please verify the correct "
                f"organism/species was selected for this dataset upload. If it was correct, "
                f"contact the gEAR team for help and reference share ID {self.share_uid}."
            )

        self._update_progress(50, "Writing initial H5AD file...")
        if adata.X is None:
            # TODO: This is currently not an option in the UI, but was suggested to be one by @jorvis
            adata = SeuratUploader.layer_to_X(adata, layer_name='data')
        h5ad_path = self.staging_area / f"{self.share_uid}.new.h5ad"
        try:
            adata.write(h5ad_path)

            # Replace the original file with the sanitized one
            #seurat_filepath.unlink()
            Path(adata_filepath).unlink()
            h5ad_path.rename(self.staging_area / f"{self.share_uid}.h5ad")
        except Exception as e:
            raise ProcessingError(
                f"Internal error while writing the initial H5AD file: {e}. Please contact the "
                f"gEAR team to resolve this issue (share ID: {self.share_uid})."
            )

        self._update_progress(65, "Seurat processing complete.")
        return h5ad_path

    def _process_excel(self) -> Path:
        """Process Excel file with expression, observations, and genes sheets."""
        self._update_progress(5, "Reading Excel file...")

        filepath = self.staging_area / f"{self.share_uid}.xlsx"

        try:
            exp_df = pd.read_excel(filepath, sheet_name='expression', index_col=0).transpose()
        except ValueError:
            raise ProcessingError(
                "No 'expression' sheet found in your Excel file. Please add a sheet named "
                "'expression' containing the expression matrix and re-upload."
            )

        try:
            X = exp_df.to_numpy()[:, 0:].astype(float)
        except ValueError:
            raise ProcessingError(
                "The 'expression' sheet contains values that aren't numbers. Please ensure every "
                "cell in the expression matrix (aside from row/column headers) is numeric."
            )

        number_obs_from_exp, number_genes_from_exp = X.shape

        self._update_progress(20, "Reading observations sheet...")

        try:
            obs_df = pd.read_excel(filepath, sheet_name='observations', index_col='observations')
        except ValueError:
            raise ProcessingError(
                "No 'observations' sheet found in your Excel file. Please add a sheet named "
                "'observations' with an 'observations' index column and re-upload."
            )

        # Validate observations
        number_obs = len(obs_df)
        if number_obs != number_obs_from_exp:
            raise ProcessingError(
                f"The 'observations' sheet has {number_obs} row(s), but the 'expression' sheet has "
                f"{number_obs_from_exp} observation row(s). Please make sure both sheets describe "
                "the same set of observations and re-upload."
            )

        if not obs_df.index.equals(exp_df.index):
            raise ProcessingError(
                "The 'observations' sheet's index doesn't match the 'expression' sheet's "
                "observation names/order. Please make sure both sheets list observations with "
                "identical names in the same order, then re-upload."
            )

        self._update_progress(35, "Reading genes sheet...")

        genes_df = self._extract_genes_sheet(filepath, exp_df)

        # Validate gene_symbol column
        if 'gene_symbol' not in genes_df.columns:
            raise ProcessingError(
                "No 'gene_symbol' column found in the genes sheet. Please add a 'gene_symbol' "
                "column listing a gene symbol for each gene and re-upload."
            )

        digit_count = genes_df['gene_symbol'].str.isnumeric().sum()
        if digit_count > 0:
            raise ProcessingError(
                f"{digit_count} value(s) in the 'gene_symbol' column are numbers rather than gene "
                "symbols (e.g. a numeric gene ID instead of a name like 'Actb'). Please provide "
                "valid gene symbols and re-upload."
            )

        # Validate gene count
        number_genes = len(genes_df)
        if number_genes != number_genes_from_exp:
            raise ProcessingError(
                f"The genes sheet has {number_genes} row(s), but the 'expression' sheet has "
                f"{number_genes_from_exp} gene column(s). Please make sure both sheets describe "
                "the same set of genes and re-upload."
            )

        if not genes_df.index.equals(exp_df.columns):
            raise ProcessingError(
                "The genes sheet's index doesn't match the 'expression' sheet's gene names/order. "
                "Please make sure both sheets list genes with identical names in the same order, "
                "then re-upload."
            )

        self._update_progress(50, "Creating AnnData object...")

        adata = anndata.AnnData(X=X, obs=obs_df, var=genes_df)

        self._update_progress(55, "Standardizing and flagging observation metadata...")
        self._sanitize_and_flag_obs_columns(adata)

        self._update_progress(60, "Writing H5AD file...")

        h5ad_path = self.staging_area / f"{self.share_uid}.h5ad"
        adata.write(h5ad_path, compression='gzip')

        self._update_progress(65, "Excel processing complete.")
        return h5ad_path

    def _extract_genes_sheet(self, filepath: Path, exp_df: pd.DataFrame) -> pd.DataFrame:
        """Extract genes sheet, falling back to expression sheet if needed."""
        try:
            genes_df = pd.read_excel(
                filepath, sheet_name='genes', index_col=0,
                converters={'gene_symbol': str}
            )
        except ValueError:
            try:
                # Try getting genes from expression sheet
                genes_df = pd.read_excel(
                    filepath, sheet_name='expression', index_col=0,
                    usecols=[0, 1]
                )
                genes_df = genes_df.drop(genes_df.columns[0], axis=1)
            except Exception:
                raise ProcessingError(
                    "No 'genes' sheet found in your Excel file, and a gene list could not be "
                    "derived from the 'expression' sheet either. Please add a sheet named 'genes' "
                    "(with at least a 'gene_symbol' column) and re-upload."
                )

        return genes_df

    def _process_mex(self) -> Path:
        """Process MEX format (matrix.mtx, barcodes.tsv, genes.tsv or the gzipped v3 equivalents)."""
        self._update_progress(5, "Reading MEX files...")

        # The files may be inside a folder within the archive
        matrix_files = sorted(self.staging_area.rglob("matrix.mtx*"))
        if not matrix_files:
            raise ProcessingError(
                "No matrix.mtx file found in your archive. Please include it and re-upload."
            )
        mex_dir = matrix_files[0].parent

        try:
            adata = sc.read_10x_mtx(mex_dir, var_names="gene_ids", cache=False)
        except Exception as e:
            raise ProcessingError(
                f"Could not read the MEX files in your archive: {e}. Please check that matrix.mtx, "
                "barcodes.tsv and genes.tsv (or matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz) "
                "are valid Cell Ranger output files and re-upload."
            )

        # gEAR expects Ensembl IDs as the var index and a 'gene_symbol' column
        adata.var = adata.var.rename(columns={"gene_symbols": "gene_symbol"})

        self._update_progress(40, "Standardizing and flagging observation metadata...")
        self._sanitize_and_flag_obs_columns(adata)

        self._update_progress(50, "Writing H5AD file...")

        h5ad_path = self.staging_area / f"{self.share_uid}.h5ad"
        adata.write(h5ad_path, compression='gzip')

        self._update_progress(65, "MEX processing complete.")
        return h5ad_path

    def _read_expression_matrix_chunks(
        self, filepath: Path, chunk_size: int, total_rows: int
    ) -> sparse.csr_matrix:
        """Read expression matrix in chunks with progress updates."""
        total_chunks = (total_rows + chunk_size - 1) // chunk_size

        last_error: Exception | None = None
        for cleanup in (False, True):
            try:
                expression_matrix = self._read_chunks_with_cleanup(
                    filepath, chunk_size, total_rows, total_chunks, cleanup=cleanup
                )
                return sparse.csr_matrix(sparse.vstack(expression_matrix))
            except Exception as e:
                last_error = e

        raise ProcessingError(
            f"Could not parse the expression matrix file: {last_error}. Please verify it is a "
            "tab-delimited file with numeric values (aside from row/column headers) and "
            f"{total_rows} data row(s) matching the genes/observations sheets, then re-upload."
        )

    def _read_chunks_with_cleanup(
        self,
        filepath: Path,
        chunk_size: int,
        total_rows: int,
        total_chunks: int,
        cleanup: bool = False,
    ) -> list:
        """
        Read expression matrix in chunks, optionally with data cleanup.

        Args:
            filepath: Path to expression matrix file
            chunk_size: Number of rows per chunk
            total_rows: Total number of rows
            total_chunks: Total number of chunks
            cleanup: Whether to apply data cleanup (whitespace, type conversion)

        Returns:
            List of sparse matrices.

        Raises:
            Whatever exception pandas/scipy raise while parsing a chunk — the caller
            (_read_expression_matrix_chunks) is responsible for retrying and reporting.
        """
        expression_matrix = []
        rows_read = 0
        reader = pd.read_csv(filepath, sep='\t', index_col=0, chunksize=chunk_size)

        for chunk_idx, chunk in enumerate(reader, 1):
            if cleanup:
                chunk = clean_chunk(chunk)

            rows_read += len(chunk)
            pct = int((rows_read / total_rows) * 100)
            self._update_progress(
                25 + round(pct * 0.4),
                f"Processing chunk {chunk_idx}/{total_chunks}"
            )
            expression_matrix.append(sparse.csr_matrix(chunk.values))

        return expression_matrix

    def _update_var_with_ensembl_ids(self, var_df: pd.DataFrame) -> pd.DataFrame:
        """Update var dataframe with Ensembl IDs, falling back to MyGene for genes the local DB misses."""
        metadata_file = self.staging_area / 'metadata.json'
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        sample_taxid = metadata.get("sample_taxid", None)
        organism_id = geardb.get_organism_id_by_taxon_id(sample_taxid)

        if not organism_id:
            raise ProcessingError(
                "Could not determine the organism for this dataset from the sample taxonomic ID "
                "in its metadata. Please verify the correct organism/species was selected for "
                "this dataset upload, then re-upload."
            )

        id_prefix = "UNMAPPED_"
        var_df = update_var_with_ensembl_ids(var_df, organism_id, id_prefix)

        unmapped_mask = var_df.index.str.startswith(id_prefix)
        if unmapped_mask.any() and sample_taxid:
            unmapped_symbols = var_df.loc[unmapped_mask, "gene_symbol"].tolist()
            try:
                mygene_mapping = map_gene_symbols_via_mygene(unmapped_symbols, sample_taxid)
            except Exception:
                mygene_mapping = {}

            if mygene_mapping:
                new_index = var_df.index.to_series()
                resolved = var_df.loc[unmapped_mask, "gene_symbol"].map(mygene_mapping)
                new_index.loc[unmapped_mask] = resolved.where(resolved.notna(), new_index.loc[unmapped_mask])
                var_df.index = pd.Index(new_index, name=var_df.index.name)

        return var_df

    def _sanitize_and_flag_obs_columns(self, adata: anndata.AnnData) -> None:
        """
        Standardize/sanitize the obs table, then scan it for numeric columns
        that look like they may actually be categorical (e.g. replicate/slide
        numbers), and record them in metadata.json for the uploader's "review
        column types" step. Operates on the in-memory AnnData object, right
        before it's written, so the H5AD never needs to be reopened from disk
        for this.
        """
        adata.obs = standardize_and_sanitize_obs(adata.obs)
        questionable = flag_ambiguous_obs_columns(adata.obs)

        metadata_file = self.staging_area / 'metadata.json'
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        metadata['questionable_obs_columns'] = questionable
        # Nothing flagged -- nothing for the user to review, so the uploader
        # can skip straight past that step.
        metadata['obs_dtype_reviewed'] = not bool(questionable)

        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=4)

    def _cancelled_result(self) -> dict:
        """Result returned when the upload was deleted while it was being processed."""
        return {
            "success": 0,
            "cancelled": True,
            "message": f"The upload was deleted (staging directory {self.staging_area} no longer exists); processing stopped.",
        }

    def _update_progress(self, progress: int, message: str) -> None:
        """Update progress and write status file; stops processing if the upload was deleted."""
        raise_if_upload_deleted(self.staging_area)
        progress = max(0, min(100, progress))  # Clamp to 0-100
        self.status['progress'] = progress
        self.status['message'] = message
        self._write_status_file()

    def _update_status(self, status_name: str, message: str) -> None:
        """Update status and write status file."""
        self.status['status'] = status_name
        self.status['message'] = message
        self._write_status_file()

    def _write_status_file(self) -> None:
        """Write current status to status.json, unless the upload (staging directory) was deleted."""
        if not self.staging_area.is_dir():
            return
        write_status(self.status_file, self.status)