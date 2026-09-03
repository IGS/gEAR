"""
cosmx_reader.py - Chunked CosMx (Nanostring) reader.

Adapted from spatialdata_io.readers.cosmx.cosmx() (spatialdata_io==0.7.1,
https://github.com/scverse/spatialdata-io).
Original Code Copyright (c) scverse / spatialdata-io contributors (BSD 3-Clause License).

MODIFICATIONS & RATIONALE:
The upstream reader loads the entire exprMat_file.csv into memory as a single dense
pandas DataFrame before converting it to a sparse matrix, which reliably OOMs on real
CosMx datasets (hundreds of thousands of cells across all FOVs). This version reads
and sparsifies that file in chunks instead, bounding peak memory to one chunk rather
than the whole file. Everything else (obs handling, per-FOV affine transforms,
image/label reading, transcripts/points) is a faithful port of the upstream logic.

DEVELOPMENT & TRANSPARENCY NOTE:
- Initial chunking refactor generated via Claude Sonnet 5.
- Reviewed, audited, and verified by @adkinsrs on the gEAR team.

NOTE: If spatialdata_io's cosmx() reader changes upstream, this may need a manual
re-sync - it does not track spatialdata_io's version.
"""

import os
import re
import tempfile
from pathlib import Path
from typing import Any, Mapping, Optional

import dask.array as da
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from anndata import AnnData
from dask_image.imread import imread
from scipy import sparse
from skimage.transform import estimate_transform
from spatialdata import SpatialData
from spatialdata._logging import logger
from spatialdata.models import Image2DModel, Labels2DModel, PointsModel, TableModel
from spatialdata.transformations.transformations import Affine, Identity
from spatialdata_io._constants._constants import CosmxKeys
from spatialdata_io.readers._utils._utils import _set_reader_metadata

DEFAULT_CHUNK_SIZE = 50_000


def _read_counts_chunked(
    counts_file: Path, chunk_size: int
) -> tuple["sparse.csr_matrix", pd.Index, pd.Index]:
    """
    Read a CosMx exprMat_file.csv in chunks, building a sparse matrix incrementally
    instead of materializing the whole file as a dense DataFrame at once.

    Returns:
        (counts_matrix, row_index, gene_columns) - counts_matrix is a sparse CSR
        matrix, row_index is the "<cell_ID>_<fov>" row index matching upstream's
        convention, gene_columns is the column (gene) index.
    """
    row_chunks = []
    index_parts = []
    gene_columns = None

    reader = pd.read_csv(
        counts_file, header=0, index_col=CosmxKeys.INSTANCE_KEY, chunksize=chunk_size
    )
    for chunk in reader:
        fov_col = chunk.pop(CosmxKeys.FOV)
        chunk_index = chunk.index.astype(str).str.cat(fov_col.astype(str).to_numpy(), sep="_")

        if gene_columns is None:
            gene_columns = chunk.columns
        elif not chunk.columns.equals(gene_columns):
            raise ValueError("Inconsistent gene columns across exprMat_file.csv chunks.")

        row_chunks.append(sparse.csr_matrix(chunk.values))
        index_parts.append(chunk_index)

    if not row_chunks:
        raise ValueError(f"No data rows found in counts file: {counts_file}.")

    counts_matrix = sparse.vstack(row_chunks, format="csr")
    row_index = pd.Index(np.concatenate([idx.to_numpy() for idx in index_parts]))
    return counts_matrix, row_index, gene_columns


def read_cosmx(
    path: "str | Path",
    dataset_id: Optional[str] = None,
    transcripts: bool = False,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    imread_kwargs: Mapping[str, Any] = {},
    image_models_kwargs: Mapping[str, Any] = {},
) -> SpatialData:
    """Read CosMx Nanostring data, chunking the counts matrix read to bound peak memory."""
    path = Path(path)

    # tries to infer dataset_id from the name of the counts file
    if dataset_id is None:
        counts_files = [f for f in os.listdir(path) if str(f).endswith(CosmxKeys.COUNTS_SUFFIX)]
        if len(counts_files) == 1:
            found = re.match(rf"(.*)_{CosmxKeys.COUNTS_SUFFIX}", counts_files[0])
            if found:
                dataset_id = found.group(1)
    if dataset_id is None:
        raise ValueError("Could not infer `dataset_id` from the name of the counts file. Please specify it manually.")

    # check for file existence
    counts_file = path / f"{dataset_id}_{CosmxKeys.COUNTS_SUFFIX}"
    if not counts_file.exists():
        raise FileNotFoundError(f"Counts file not found: {counts_file}.")
    if transcripts:
        transcripts_file = path / f"{dataset_id}_{CosmxKeys.TRANSCRIPTS_SUFFIX}"
        if not transcripts_file.exists():
            raise FileNotFoundError(f"Transcripts file not found: {transcripts_file}.")
    else:
        transcripts_file = None
    meta_file = path / f"{dataset_id}_{CosmxKeys.METADATA_SUFFIX}"
    if not meta_file.exists():
        raise FileNotFoundError(f"Metadata file not found: {meta_file}.")
    fov_file = path / f"{dataset_id}_{CosmxKeys.FOV_SUFFIX}"
    if not fov_file.exists():
        raise FileNotFoundError(f"Found field of view file: {fov_file}.")
    images_dir = path / CosmxKeys.IMAGES_DIR
    if not images_dir.exists():
        raise FileNotFoundError(f"Images directory not found: {images_dir}.")
    labels_dir = path / CosmxKeys.LABELS_DIR
    if not labels_dir.exists():
        raise FileNotFoundError(f"Labels directory not found: {labels_dir}.")

    # Chunked, incrementally-sparse read of the counts matrix - this is the actual fix.
    counts_matrix, counts_index, gene_columns = _read_counts_chunked(counts_file, chunk_size)

    obs = pd.read_csv(meta_file, header=0, index_col=CosmxKeys.INSTANCE_KEY)
    obs[CosmxKeys.FOV] = pd.Categorical(obs[CosmxKeys.FOV].astype(str))
    obs[CosmxKeys.REGION_KEY] = pd.Categorical(obs[CosmxKeys.FOV].astype(str).apply(lambda s: s + "_labels"))
    obs[CosmxKeys.INSTANCE_KEY] = obs.index.astype(np.int64)
    obs.rename_axis(None, inplace=True)
    obs.index = obs.index.astype(str).str.cat(obs[CosmxKeys.FOV].to_numpy(), sep="_")

    common_index = obs.index.intersection(counts_index)

    # Map common_index -> integer row positions in the sparse counts matrix, since
    # counts_matrix is no longer a DataFrame we can .loc[] into directly.
    counts_row_positions = pd.Series(np.arange(len(counts_index)), index=counts_index)
    row_positions = counts_row_positions.loc[common_index].to_numpy()

    adata = AnnData(
        counts_matrix[row_positions, :],
        obs=obs.loc[common_index, :],
    )
    adata.var_names = gene_columns

    # Drop the redundant "cell_id" column if present - it conflicts with the
    # "cell_ID" instance key and otherwise fails TableModel validation below.
    if "cell_id" in adata.obs.columns:
        adata.obs = adata.obs.drop(columns=["cell_id"])

    table = TableModel.parse(
        adata,
        region=list(set(adata.obs[CosmxKeys.REGION_KEY].astype(str).tolist())),
        region_key=CosmxKeys.REGION_KEY.value,
        instance_key=CosmxKeys.INSTANCE_KEY.value,
    )

    fovs_counts = list(map(str, adata.obs.fov.astype(int).unique()))

    affine_transforms_to_global = {}

    for fov in fovs_counts:
        idx = table.obs.fov.astype(str) == fov
        loc = table[idx, :].obs[[CosmxKeys.X_LOCAL_CELL, CosmxKeys.Y_LOCAL_CELL]].values
        glob = table[idx, :].obs[[CosmxKeys.X_GLOBAL_CELL, CosmxKeys.Y_GLOBAL_CELL]].values
        out = estimate_transform(ttype="affine", src=loc, dst=glob)
        affine_transforms_to_global[fov] = Affine(
            out.params,
            input_axes=("x", "y"),
            output_axes=("x", "y"),
        )

    table.obsm["global"] = table.obs[[CosmxKeys.X_GLOBAL_CELL, CosmxKeys.Y_GLOBAL_CELL]].to_numpy()
    table.obsm["spatial"] = table.obs[[CosmxKeys.X_LOCAL_CELL, CosmxKeys.Y_LOCAL_CELL]].to_numpy()
    table.obs.drop(
        columns=[CosmxKeys.X_LOCAL_CELL, CosmxKeys.Y_LOCAL_CELL, CosmxKeys.X_GLOBAL_CELL, CosmxKeys.Y_GLOBAL_CELL],
        inplace=True,
    )

    # prepare to read images and labels
    file_extensions = (".jpg", ".png", ".jpeg", ".tif", ".tiff")
    pat = re.compile(r".*_F(\d+)")

    # check if fovs are correct for images and labels
    fovs_images = []
    for fname in os.listdir(path / CosmxKeys.IMAGES_DIR):
        if fname.endswith(file_extensions):
            fovs_images.append(str(int(pat.findall(fname)[0])))

    fovs_labels = []
    for fname in os.listdir(path / CosmxKeys.LABELS_DIR):
        if fname.endswith(file_extensions):
            fovs_labels.append(str(int(pat.findall(fname)[0])))

    fovs_images_and_labels = set(fovs_images).intersection(set(fovs_labels))
    fovs_diff = fovs_images_and_labels.difference(set(fovs_counts))
    if len(fovs_diff):
        logger.warning(
            f"Found images and labels for {len(fovs_images)} FOVs, but only {len(fovs_counts)} FOVs in the counts file.\n"
            + f"The following FOVs are missing: {fovs_diff} \n"
            + "... will use only fovs in Table."
        )

    # read images
    images = {}
    for fname in os.listdir(path / CosmxKeys.IMAGES_DIR):
        if fname.endswith(file_extensions):
            fov = str(int(pat.findall(fname)[0]))
            if fov in fovs_counts:
                aff = affine_transforms_to_global[fov]
                im = imread(path / CosmxKeys.IMAGES_DIR / fname, **imread_kwargs).squeeze()
                flipped_im = da.flip(im, axis=0)
                parsed_im = Image2DModel.parse(
                    flipped_im,
                    transformations={
                        fov: Identity(),
                        "global": aff,
                        "global_only_image": aff,
                    },
                    dims=("y", "x", "c"),
                    rgb=None,
                    **image_models_kwargs,
                )
                images[f"{fov}_image"] = parsed_im
            else:
                logger.warning(f"FOV {fov} not found in counts file. Skipping image {fname}.")

    # read labels
    labels = {}
    for fname in os.listdir(path / CosmxKeys.LABELS_DIR):
        if fname.endswith(file_extensions):
            fov = str(int(pat.findall(fname)[0]))
            if fov in fovs_counts:
                aff = affine_transforms_to_global[fov]
                la = imread(path / CosmxKeys.LABELS_DIR / fname, **imread_kwargs).squeeze()
                flipped_la = da.flip(la, axis=0)
                parsed_la = Labels2DModel.parse(
                    flipped_la,
                    transformations={
                        fov: Identity(),
                        "global": aff,
                        "global_only_labels": aff,
                    },
                    dims=("y", "x"),
                    **image_models_kwargs,
                )
                labels[f"{fov}_labels"] = parsed_la
            else:
                logger.warning(f"FOV {fov} not found in counts file. Skipping labels {fname}.")

    points: dict = {}
    if transcripts:
        with tempfile.TemporaryDirectory() as tmpdir:
            assert transcripts_file is not None
            transcripts_data = pd.read_csv(transcripts_file, header=0)
            transcripts_data.to_parquet(Path(tmpdir) / "transcripts.parquet")

            ptable = pq.read_table(Path(tmpdir) / "transcripts.parquet")
            for fov in fovs_counts:
                aff = affine_transforms_to_global[fov]
                sub_table = ptable.filter(pa.compute.equal(ptable.column(CosmxKeys.FOV), int(fov))).to_pandas()
                sub_table[CosmxKeys.INSTANCE_KEY] = sub_table[CosmxKeys.INSTANCE_KEY].astype("category")
                # we rename z because we want to treat the data as 2d
                sub_table.rename(columns={"z": "z_raw"}, inplace=True)
                if len(sub_table) > 0:
                    points[f"{fov}_points"] = PointsModel.parse(
                        sub_table,
                        coordinates={"x": CosmxKeys.X_LOCAL_TRANSCRIPT, "y": CosmxKeys.Y_LOCAL_TRANSCRIPT},
                        feature_key=CosmxKeys.TARGET_OF_TRANSCRIPT,
                        instance_key=CosmxKeys.INSTANCE_KEY,
                        transformations={
                            fov: Identity(),
                            "global": aff,
                            "global_only_labels": aff,
                        },
                    )

    sdata = SpatialData(images=images, labels=labels, points=points, tables={"table": table})
    return _set_reader_metadata(sdata, "cosmx")
