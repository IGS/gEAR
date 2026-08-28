# utils.py - Functions that are used across multiple scripts.

# Some of these originally started out as code from individual scripts, such as the "bin" directory,
# but were moved to a common location to be shared across multiple scripts.

import functools
import os
import sys
import typing

if typing.TYPE_CHECKING:
    import pandas as pd
    from anndata import AnnData

def set_memory_limit_from_cgroup(fraction: float = 0.9) -> None:
    """
    Self-impose a process memory ceiling based on the container's cgroup memory limit,
    so that approaching the real limit raises a catchable MemoryError instead of the
    kernel OOM-killer sending an uncatchable SIGKILL (exit code 137).

    Reads the cgroup v2 limit first (/sys/fs/cgroup/memory.max), falling back to
    cgroup v1 (/sys/fs/cgroup/memory/memory.limit_in_bytes). If no bounded limit is
    found (unbounded, or the files aren't present), this is a no-op and processes
    remain subject to the OS OOM-killer as before.

    This is best-effort defensive setup: any failure here is caught and logged
    rather than raised, so it can never prevent a caller from starting up.

    Parameters
    ----------
    fraction : float, optional (default: 0.9)
        Fraction of the detected cgroup memory limit to set as the RLIMIT_AS ceiling.
    """
    import resource

    # Cgroup v1 reports an implausibly large number (close to 2**63, platform max)
    # to mean "no limit" rather than a sentinel string like v2's "max".
    UNBOUNDED_V1_THRESHOLD = 2**62

    limit_bytes = None

    try:
        cgroup_v2_path = "/sys/fs/cgroup/memory.max"
        cgroup_v1_path = "/sys/fs/cgroup/memory/memory.limit_in_bytes"

        if os.path.exists(cgroup_v2_path):
            with open(cgroup_v2_path) as f:
                value = f.read().strip()
            if value != "max":
                limit_bytes = int(value)
        elif os.path.exists(cgroup_v1_path):
            with open(cgroup_v1_path) as f:
                value = int(f.read().strip())
            if value < UNBOUNDED_V1_THRESHOLD:
                limit_bytes = value

        if limit_bytes is None:
            print(
                "set_memory_limit_from_cgroup: no bounded cgroup memory limit found; "
                "not setting a self-imposed RLIMIT_AS ceiling.",
                file=sys.stderr,
            )
            return

        target = int(limit_bytes * fraction)
        resource.setrlimit(resource.RLIMIT_AS, (target, target))
        print(
            f"set_memory_limit_from_cgroup: cgroup limit is {limit_bytes} bytes; "
            f"set RLIMIT_AS to {target} bytes ({fraction:.0%}).",
            file=sys.stderr,
        )
    except Exception as e:
        print(f"set_memory_limit_from_cgroup: failed to set memory limit: {e}", file=sys.stderr)


def catch_memory_error() -> typing.Callable:
    """
    A decorator factory that catches MemoryError exceptions in the decorated function.

    Returns:
        Callable: A decorator that wraps the target function. If a MemoryError is raised during
        execution, it prints an error message to stderr and returns a tuple containing a result
        dictionary and a 500 status code.

    Example:
        @catch_memory_error()
        def my_function():
            # function implementation
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except MemoryError as e:
                print(f"Exceeded memory in {func.__name__}: {e}", file=sys.stderr)

                result = {
                    "message": "Exceeded memory limit",
                    "success": -1,
                    "error": str(e),
                }

                return result, 500

        return wrapper

    return decorator


def update_adata_with_ensembl_ids(
    adata: "AnnData", organism: int, id_prefix: str, verbose: bool = False
) -> "AnnData":
    """
    Updates the gene identifiers in an AnnData object to Ensembl IDs by mapping gene symbols to Ensembl IDs
    using a database lookup for the specified organism and a set of Ensembl releases. The function selects
    the Ensembl release with the highest number of matches to the gene symbols in the input AnnData object.

    Genes that cannot be mapped to an Ensembl ID are retained with their original identifiers, prefixed by
    the provided `id_prefix`. The function ensures that the shape and metadata of the AnnData object are
    preserved, and handles duplicate gene symbols by dropping them prior to mapping.

    Parameters
    ----------
    adata : AnnData
        The input AnnData object containing gene expression data. Gene symbols are expected to be in `adata.var.index`.
    organism : int
        The organism primary key ID in the `geardb` database as found in the "organism" table.
    id_prefix : str
        Prefix to use for genes that cannot be mapped to an Ensembl ID.
    verbose : bool, optional (default: False)
        If True, prints detailed progress and debugging information.

    Returns
    -------
    AnnData
        A new AnnData object with gene identifiers updated to Ensembl IDs where possible. Unmapped genes are
        retained with their original identifiers, prefixed by `id_prefix`. All original metadata and structure
        are preserved.

    Notes
    -----
    - Requires access to a `geardb` database with gene symbol to Ensembl ID mappings.
    - Drops duplicate gene symbols prior to mapping to ensure one-to-one mapping.
    - Preserves all AnnData fields (obs, obsm, obsp, varm, varp, uns) in the output.
    """

    import anndata as ad
    import geardb
    import pandas as pd

    (_, n_genes) = adata.shape
    ensembl_releases = [84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94]

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    query = """
        SELECT ensembl_id, gene_symbol
          FROM gene
         WHERE organism_id = %s
           AND ensembl_release = %s
    """

    best_release = None
    best_count = 0
    best_df = None

    # There are some cases where there are duplicate gene symbol,
    # we do not want to keep these so we drop them here in order
    # preserve obs shape when joining ensembl ids and resassigning
    # var later on.
    duplicated_genes = adata.var.index.duplicated()

    if verbose:
        print(f"Duplicated Genes: {duplicated_genes.sum()}")

    # print(f"The file, {args.input_file} has {duplicated_genes.sum()} duplicate genes. These will be dropped.")
    adata = adata[:, ~duplicated_genes]

    if verbose:
        print("\nOriginal loaded adata\n")
        print(adata)

    # If adata.var is an empty dataframe, make note that the index is the original gene symbol column
    # Ensures the `adata_unmapped_var` rename aligns with the original gene symbol column in adata.var
    orig_gene_column = "genes"
    if adata.var.empty:
        orig_gene_column = "index"

    for release in ensembl_releases:
        if verbose:
            print(
                "INFO: comparing with ensembl release: {0} ... ".format(release), end=""
            )
        cursor.execute(query, (organism, release))

        df = pd.DataFrame(cursor.fetchall(), columns=cursor.column_names)

        # Query can return different ensembl ids for a gene symbol,
        # we want to drop duplicates so we have only one ensembl id
        # per gene symbol
        df = df.drop_duplicates(subset=["gene_symbol"])
        df = df.set_index("gene_symbol")

        merged_df = adata.var.join(df, how="inner")
        (row_count, _) = merged_df.shape

        if verbose:
            print(" found {0} matches".format(row_count))

        if row_count > best_count:
            best_count = row_count
            best_release = release
            best_df = merged_df

    if best_df is None:
        raise ValueError(
            "No matches found for any of the Ensembl releases. Please check your input data and organism ID."
        )

    if verbose:
        print(f"\nBest release: {best_release}")
        print(f"Matches for release: {best_count}")
        print(f"Original # Genes: {n_genes}")
        print(f"Genes lost: {n_genes - best_count}\n")

    # Now we have our best release and ensembl ids for those gene symbols,

    # Get separate adata for those where the gene symbols were mapped and where they weren't
    genes_present_filter = adata.var.index.isin(best_df.index)
    adata_present = adata[:, genes_present_filter]
    adata_not_present = adata[:, ~genes_present_filter]

    # If the data already had a 'gene symbol' let's rename it
    if "gene_symbol" in best_df.columns:
        if verbose:
            print(
                "WARN: Found gene_symbol column already in dataset, renaming to gene_symbol_original"
            )
        best_df = best_df.rename(columns={"gene_symbol": "gene_symbol_original"})

    ensembl_id_var = (
        best_df.reset_index()
        .rename(columns={"index": "gene_symbol"})
        .set_index("ensembl_id")
    )

    if verbose:
        print("ENSEMBL_ID_VAR")
        print(ensembl_id_var)

    # Currently creating a new AnnData object because of
    # trouble getting adata.var = merged_var to persist
    adata_with_ensembl_ids = ad.AnnData(
        adata_present.X,
        obs=adata_present.obs,
        var=ensembl_id_var,
        # May not use these directly, but need to pass them through to preserve them
        # Note: there are other fields, like layers, that could be added here if needed
        obsm=adata_present.obsm,
        obsp=adata_present.obsp,
        varm=adata_present.varm,
        varp=adata_present.varp,
        uns=adata_present.uns,
    )

    if verbose:
        print(adata_with_ensembl_ids.obs.columns)

    ## Now combine the unmapped dataframe with this one, first making the needed edits
    if "gene_symbol" in adata_not_present.var.columns:
        adata_not_present.var = adata_not_present.var.rename(
            columns={"gene_symbol": "gene_symbol_original"}
        )

    if verbose:
        print("ADATA_NOT_PRESENT.VAR")
        print(adata_not_present.var)
        print(adata_not_present.obs.columns)

    # Splitting code over multiple lines requires a "\" at the end.
    adata_unmapped_var = (
        adata_not_present.var.reset_index(names=orig_gene_column)
        .rename(columns={orig_gene_column: "gene_symbol"})
        .set_index(id_prefix + adata_not_present.var.index.astype(str))
    )

    adata_unmapped = ad.AnnData(
        X=adata_not_present.X,
        obs=adata_with_ensembl_ids.obs,
        var=adata_unmapped_var,
        obsm=adata_not_present.obsm,
        obsp=adata_not_present.obsp,
        varm=adata_not_present.varm,
        varp=adata_not_present.varp,
        uns=adata_not_present.uns,
    )
    adata_unmapped.var.index.name = "ensembl_id"

    if verbose:
        print("ADATA UNMAPPED.VAR")
        print(adata_unmapped.var)

    # Concatenate the two AnnData objects (adata_with_ensembl_ids and adata_unmapped)
    # This will leave an index-only set of observations (WHY?) so we need to reassign the obs.
    # Reassign the other structures as well to ensure they are preserved.

    adata = ad.concat([adata_with_ensembl_ids, adata_unmapped], axis=1)
    adata.obs = adata_present.obs
    adata.obsm = adata_present.obsm
    adata.obsp = adata_present.obsp
    adata.varp = adata_present.varp
    adata.uns = adata_present.uns

    if verbose:
        print("ADATA CONCAT")
        print(adata)
        # print("VAR CONCAT")
        # print(adata.var)
        # print("OBS_CONCAT")
        # print(adata.obs)
        # print(adata.X)

    return adata

def update_var_with_ensembl_ids(
    var_df: "pd.DataFrame", organism: int, id_prefix: str, verbose: bool = False
) -> "pd.DataFrame":
    """
    Updates gene identifiers in a var dataframe to Ensembl IDs.

    Designed to work with backed-mode AnnData objects where only metadata
    needs to be updated without loading the full expression matrix.

    Parameters
    ----------
    var_df : pd.DataFrame
        The var dataframe from an AnnData object (gene-level metadata).
    organism : int
        Organism primary key ID in geardb database.
    id_prefix : str
        Prefix for unmapped genes.
    verbose : bool, optional (default: False)
        Print progress information.

    Returns
    -------
    pd.DataFrame
        Updated var dataframe with Ensembl IDs as index.
    """
    import geardb
    import pandas as pd

    ensembl_releases = [84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94]

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    query = """
        SELECT ensembl_id, gene_symbol
          FROM gene
         WHERE organism_id = %s
           AND ensembl_release = %s
    """

    # We want to ensure that the var_df does not contain duplicate gene symbols, as this would lead to ambiguity in mapping.
    # Not only that, dropping duplicates would cause a mismatch with the shape of Anndata.X
    if var_df.index.duplicated().any():
        raise ValueError(
            "var_df has duplicate gene symbols; drop duplicates (and the matching "
            "X columns) before calling update_var_with_ensembl_ids."
        )

    # Preserve the caller's original row order through the split/concat below,
    # since AnnData.var assignment is positional and must match X's column order.
    var_df = var_df.copy()
    var_df["_orig_order"] = range(len(var_df))
    orig_gene_column = "index"

    best_release = None
    best_count = 0
    best_df = None

    for release in ensembl_releases:
        if verbose:
            print(f"Comparing with Ensembl release {release} ... ", end="")

        cursor.execute(query, (organism, release))
        df = pd.DataFrame(cursor.fetchall(), columns=cursor.column_names)
        df = df.drop_duplicates(subset=["gene_symbol"])
        df = df.set_index("gene_symbol")

        merged_df = var_df.join(df, how="inner")
        row_count = len(merged_df)

        if verbose:
            print(f"found {row_count} matches")

        if row_count > best_count:
            best_count = row_count
            best_release = release
            best_df = merged_df

    if best_df is None:
        raise ValueError("No Ensembl matches found for organism.")

    if verbose:
        print(f"Best release: {best_release}, Matches: {best_count}")

    # Mapped genes
    genes_present_filter = var_df.index.isin(best_df.index)
    var_present = var_df[genes_present_filter]
    var_not_present = var_df[~genes_present_filter]

    # Rename existing gene_symbol column if present
    if "gene_symbol" in best_df.columns:
        best_df = best_df.rename(columns={"gene_symbol": "gene_symbol_original"})

    # Update mapped genes
    ensembl_id_var = (
        best_df.reset_index()
        .rename(columns={"index": "gene_symbol"})
        .set_index("ensembl_id")
    )
    ensembl_id_var.index.name = "ensembl_id"

    # Update unmapped genes
    if "gene_symbol" in var_not_present.columns:
        var_not_present = var_not_present.rename(
            columns={"gene_symbol": "gene_symbol_original"}
        )

    unmapped_var = (
        var_not_present.reset_index(names=orig_gene_column)
        .rename(columns={orig_gene_column: "gene_symbol"})
        .set_index(id_prefix + var_not_present.index.astype(str))
    )
    unmapped_var.index.name = "ensembl_id"

    # Combine and restore order
    result = pd.concat([ensembl_id_var, unmapped_var], axis=0)
    result = result.sort_values("_orig_order").drop(columns="_orig_order")

    return result

def map_gene_symbols_via_mygene(
    gene_symbols: list, taxid: int, verbose: bool = False
) -> dict:
    """
    Query MyGene.info for Ensembl gene IDs matching the given gene symbols.

    Intended as a fallback for genes the local `gene` table can't resolve
    (e.g. older Ensembl releases). Best-effort: symbols MyGene can't resolve,
    or any request failure after retries, should be handled by the caller.

    Parameters
    ----------
    gene_symbols : list
        Gene symbols to look up.
    taxid : int
        NCBI taxonomy ID for the organism.
    verbose : bool, optional (default: False)
        Print progress information.

    Returns
    -------
    dict
        Mapping of gene_symbol -> ensembl_id for symbols MyGene resolved.
        Unresolved symbols are omitted.
    """
    import time

    import mygene

    if not gene_symbols:
        return {}

    max_retries = 5
    base_delay = 2  # seconds

    mg_genes = None
    for attempt in range(max_retries):
        try:
            mg = mygene.MyGeneInfo()
            mg_genes = mg.querymany(gene_symbols, scopes="symbol", fields="ensembl.gene", species=str(taxid))
            break
        except Exception as e:
            is_server_error = "500" in str(e) or "internal server error" in str(e).lower()
            if is_server_error and attempt < max_retries - 1:
                delay = base_delay ** (attempt + 1)
                if verbose:
                    print(
                        f"MyGene API returned a 500 error (attempt {attempt + 1}/{max_retries}). "
                        f"Retrying in {delay}s..."
                    )
                time.sleep(delay)
            else:
                if verbose:
                    print(f"Error occurred while querying MyGene: {e}")
                raise

    mapping = {}
    for mg_gene in mg_genes or []:
        if "ensembl" not in mg_gene:
            continue
        ensembl = mg_gene["ensembl"]
        mapping[mg_gene["query"]] = ensembl[0]["gene"] if isinstance(ensembl, list) else ensembl["gene"]

    return mapping
