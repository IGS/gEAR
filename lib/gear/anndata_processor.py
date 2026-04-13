"""
AnndataProcessor - Process expression datasets in various formats to H5AD.

Handles conversion, sanitization, and primary analysis of expression datasets
in formats: H5AD, 3-tab, Excel, MEX.
"""

import gc
import json
import os
from pathlib import Path
import tarfile
import zipfile

import anndata
import pandas as pd
from scipy import sparse

import geardb
from gear.primary_analysis import add_primary_analysis_to_dataset, PrimaryAnalysisProcessingError
from gear.utils import update_var_with_ensembl_ids



def write_status(status_file, status):
    """Write status dictionary to JSON file."""
    with open(status_file, 'w') as f:
        f.write(json.dumps(status, indent=2))

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

def sanitize_obs_for_h5ad(obs_df: pd.DataFrame) -> pd.DataFrame:
    """Sanitize observation dataframe for H5AD storage."""
    for col in obs_df.columns:
        if obs_df[col].dtype == 'object':
            obs_df[col] = obs_df[col].fillna('').astype(str)
    return obs_df

def categorize_observation_columns(obs: pd.DataFrame) -> None:
    """Categorize and convert specific observation columns."""
    for str_type in ['cell_type', 'condition', 'time_point', 'time_unit']:
        if str_type in obs.columns:
            obs[str_type] = pd.Categorical(obs[str_type])

    for num_type in ['replicate', 'time_point_order']:
        if num_type in obs.columns:
            obs[num_type] = pd.to_numeric(obs[num_type])


def package_content_type(filenames: list[str]) -> str | None:
        #print("DEBUG: filenames", file=sys.stderr, flush=True)
        #print(filenames, file=sys.stderr, flush=True)
        """
        mex:
        matrix.mtx
        barcodes.tsv
        genes.tsv

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

        if 'matrix.mtx' in filenames and 'barcodes.tsv' in filenames and 'genes.tsv' in filenames:
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
        self.job_id = job_id
        self.share_uid = share_uid
        self.staging_area = staging_area
        self.status_file = status_file
        self.dataset_uid = dataset_uid
        self.status = {
            "process_id": os.getpid(),
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
            h5ad_path = self._process_by_format(dataset_format)

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
                    raise ProcessingError(f"Primary analysis failed: {str(e)}")

            self._update_progress(100, "Dataset processed successfully.")
            return {"success": 1, "message": "Dataset processed successfully."}

        except ProcessingError as e:
            self._update_status("error", str(e))
            return {"success": 0, "message": str(e)}
        except Exception as e:
            import traceback
            traceback.print_exc()
            self._update_status("error", f"Unexpected error: {str(e)}")
            return {"success": 0, "message": f"Unexpected error: {str(e)}"}
        finally:
            gc.collect()

    def _process_by_format(self, dataset_format: str) -> Path:
        """Route to appropriate processor based on dataset format."""
        if dataset_format == "h5ad":
            return self._process_h5ad()
        elif dataset_format == "mex_3tab":
            return self._process_mex_3tab()
        elif dataset_format == "excel":
            return self._process_excel()
        elif dataset_format == "rdata":
            raise NotImplementedError("RData format processing not yet implemented.")
            return self._process_rdata()
        else:
            raise ProcessingError(f"Unsupported dataset format: {dataset_format}")

    def _process_rdata(self):
        pass

    def _process_h5ad(self) -> Path:
        """Process .h5ad file with backed mode for memory efficiency."""
        self._update_progress(5, "Reading H5AD file...")

        filepath = self.staging_area / f"{self.share_uid}.h5ad"

        # Read in backed mode to avoid loading full file into memory
        adata = anndata.read_h5ad(filepath, backed='r')

        self._update_progress(15, "Sanitizing observation metadata...")

        # Convert obs to memory (typically small)
        obs = adata.obs
        categorize_observation_columns(obs)
        obs = sanitize_obs_for_h5ad(obs)

        # Update var if gene_symbol missing
        var = adata.var if "gene_symbol" not in adata.var.columns else None
        if var is not None:
            self._update_progress(25, "Mapping gene symbols via Ensembl...")
            var = self._update_var_with_ensembl_ids(var)

        # Close backed file before final write
        adata.file.close()

        self._update_progress(50, "Writing H5AD file...")

        # Re-read for full write (necessary for consistency)
        adata = anndata.read_h5ad(filepath)
        adata.obs = obs
        if var is not None:
            adata.var = var

        h5ad_temp = self.staging_area / f"{self.share_uid}.new.h5ad"
        adata.write(h5ad_temp, compression='gzip')

        filepath.unlink()
        h5ad_temp.rename(filepath)

        self._update_progress(65, "H5AD processing complete.")
        return filepath

    def _process_mex_3tab(self) -> Path:
        # Extract the file
        compression_format = None
        filename = self.staging_area / f"{self.share_uid}.tar.gz"

        if filename.exists():
            compression_format = 'tarball'
        else:
            filename = self.staging_area / f"{self.share_uid}.zip"

            if filename.exists():
                compression_format = 'zip'
            else:
                raise ProcessingError("No tarball or zip file found for MEX/3-tab dataset.")

        files_extracted = []

        if compression_format == 'tarball':
            try:
                with tarfile.open(filename) as tf:
                    for entry in tf:
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
            except tarfile.ReadError:
                raise ProcessingError("Bad tarball file. Couldn't extract the tarball.")

        if compression_format == 'zip':
            try:
                with zipfile.ZipFile(filename) as zf:
                    for entry in zf.infolist():
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
            except zipfile.BadZipFile:
                raise ProcessingError("Bad zip file. Couldn't extract the zip file.")

        # Determine the dataset type
        dataset_type = package_content_type(files_extracted)

        if dataset_type is None:
            raise ProcessingError("Unsupported dataset format. Couldn't tell type from file names within the archive.")

        # Call the appropriate function
        if dataset_type == 'threetab':
            return self._process_threetab()
        elif dataset_type == 'mex':
            return self._process_mex()

        # If no valid dataset type is found, raise an error
        raise ProcessingError("Failed to process MEX/3-tab dataset. No valid dataset type identified.")


    def _process_threetab(self) -> Path:
        """Process 3-tab format (expression, genes, observations) with chunking."""
        self._update_progress(5, "Extracting and validating files...")

        # Extract/find the three required files
        expression_matrix_path, obs, var = self._extract_threetab_files()

        self._update_progress(15, "Categorizing observations...")
        categorize_observation_columns(obs)

        self._update_progress(25, "Processing expression matrix in chunks...")

        # Process expression matrix in chunks
        chunk_size = 500
        expression_matrix = self._read_expression_matrix_chunks(
            expression_matrix_path, chunk_size, total_rows=len(var)
        )

        # Create AnnData object
        adata = anndata.AnnData(X=expression_matrix, obs=var, var=obs)
        adata = adata.transpose()
        adata.obs = sanitize_obs_for_h5ad(adata.obs)

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
            raise ProcessingError("No observations file found (expected observations.tab or COLmeta.tab).")
        if var is None:
            raise ProcessingError("No genes file found (expected genes.tab or ROWmeta.tab).")
        if expression_matrix_path is None:
            raise ProcessingError("No expression file found (expected expression.tab or DataMTX.tab).")

        return expression_matrix_path, obs, var

    def _process_excel(self) -> Path:
        """Process Excel file with expression, observations, and genes sheets."""
        self._update_progress(5, "Reading Excel file...")

        filepath = self.staging_area / f"{self.share_uid}.xlsx"
        exp_df = pd.read_excel(filepath, sheet_name='expression', index_col=0).transpose()

        try:
            X = exp_df.to_numpy()[:, 0:].astype(float)
        except ValueError:
            raise ProcessingError("Encountered unexpected value type. Expected float type in expression matrix.")

        number_obs_from_exp, number_genes_from_exp = X.shape

        self._update_progress(20, "Reading observations sheet...")

        try:
            obs_df = pd.read_excel(filepath, sheet_name='observations', index_col='observations')
        except ValueError:
            raise ProcessingError("No observations sheet found. Expected spreadsheet sheet named 'observations'.")

        # Validate observations
        number_obs = len(obs_df)
        if number_obs != number_obs_from_exp:
            raise ProcessingError(
                f"Observations sheet error. Row count ({number_obs}) must match "
                f"expression sheet row count ({number_obs_from_exp})."
            )

        if not obs_df.index.equals(exp_df.index):
            raise ProcessingError(
                "Observations sheet error. Index names and order must match expression sheet rows."
            )

        self._update_progress(35, "Reading genes sheet...")

        genes_df = self._extract_genes_sheet(filepath, exp_df)

        # Validate gene_symbol column
        if 'gene_symbol' not in genes_df.columns:
            raise ProcessingError("Failed to find gene_symbol column in genes tab.")

        digit_count = genes_df['gene_symbol'].str.isnumeric().sum()
        if digit_count > 0:
            raise ProcessingError(f"{digit_count} gene symbols are numbers, not gene symbols.")

        categorize_observation_columns(obs_df)

        # Validate gene count
        number_genes = len(genes_df)
        if number_genes != number_genes_from_exp:
            raise ProcessingError(
                f"Genes sheet error. Row count ({number_genes}) must match "
                f"expression sheet row count ({number_genes_from_exp})."
            )

        if not genes_df.index.equals(exp_df.columns):
            raise ProcessingError("Genes sheet error. Index names and order must match expression sheet columns.")

        self._update_progress(50, "Creating AnnData object...")

        adata = anndata.AnnData(X=X, obs=obs_df, var=genes_df)
        adata.obs = sanitize_obs_for_h5ad(adata.obs)

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
            except Exception as err:
                raise ProcessingError(f"No 'genes' sheet found. {str(err)}")

        return genes_df

    def _process_mex(self) -> Path:
        """Process MEX format (matrix.mtx, barcodes.tsv, genes.tsv)."""
        raise NotImplementedError("MEX format processing not yet implemented.")

    def _read_expression_matrix_chunks(
        self, filepath: Path, chunk_size: int, total_rows: int
    ) -> sparse.csr_matrix:
        """Read expression matrix in chunks with progress updates."""
        total_chunks = (total_rows + chunk_size - 1) // chunk_size

        # Try reading without cleanup first
        expression_matrix = self._read_chunks_with_cleanup(
            filepath, chunk_size, total_rows, total_chunks, cleanup=False
        )

        if expression_matrix:
            return sparse.csr_matrix(sparse.vstack(expression_matrix))

        # Retry with data cleanup if first attempt fails
        expression_matrix = self._read_chunks_with_cleanup(
            filepath, chunk_size, total_rows, total_chunks, cleanup=True
        )
        if expression_matrix:
            return sparse.csr_matrix(sparse.vstack(expression_matrix))

        raise ProcessingError("Failed to read expression matrix in chunks, even after cleanup.")

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
            List of sparse matrices, empty list if reading fails
        """
        expression_matrix = []
        rows_read = 0
        reader = pd.read_csv(filepath, sep='\t', index_col=0, chunksize=chunk_size)

        try:
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

        except Exception:
            return []

    def _update_var_with_ensembl_ids(self, var_df: pd.DataFrame) -> pd.DataFrame:
        """Update var dataframe with Ensembl IDs."""
        metadata_file = self.staging_area / 'metadata.json'
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)

        sample_taxid = metadata.get("sample_taxid", None)
        organism_id = geardb.get_organism_id_by_taxon_id(sample_taxid)

        if not organism_id:
            raise ProcessingError("Could not determine organism ID from sample taxonomic ID.")

        return update_var_with_ensembl_ids(var_df, organism_id, "UNMAPPED_")

    def _update_progress(self, progress: int, message: str) -> None:
        """Update progress and write status file."""
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
        """Write current status to status.json."""
        write_status(self.status_file, self.status)