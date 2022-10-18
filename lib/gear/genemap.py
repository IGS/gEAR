import os
import pandas as pd
import sys

class GeneMap:
    """
    There are a few different use cases for mapping genes between organisms.

       - orthology mappings
       - host/pathogen interaction mappings
       - others I can't think of

    For whatever the reason, this module works based on the existence of mapping
    H5AD files, each containing a dataframe with this structure:

      - Source identifier (index)
      - Source gene symbol
      - Target identifier
      - Target gene symbol

    The inputs to be mapped can take a few forms:

      - GeneCart object
      - Dataset object
      - AnnData object
      - Single gene symbol string

    In each case, the source organism and dest organism IDs and annotation sources
    must also be defined.

    Example usage:

       gm = GeneMap(source_org_id=1, dest_org_id=1, source_annot_type='ensembl', dest_annot_type='ensembl')

       # GeneCart support
       gc = GeneCart(...)
       filtered_gc = gm.map(gc)

       # Dataset support
       ds = Dataset(...)
       filtered_ds = gm.map(ds)

       # AnnData object support
       ad = AnnData(...)
       filtered_ad = gm.map(ad)

       # Simple string support
       mapped_gene = gm.map("Pou4f3")

    Any unmapped genes are simply excluded from the returned data structure, or "None" is
    returned in the case of an unmapped string.

    """
    def __init__(self, source_org_id=None, source_annot_type=None,
                 dest_org_id=None, dest_annot_type=None):

        if source_org_id:
            self.source_org_id = source_org_id
        else:
            raise Exception("source_org_id is a required argument when creating a GeneMap")

        if source_annot_type:
            self.source_annot_type = source_annot_type
        else:
            raise Exception("source_annot_type is a required argument when creating a GeneMap")

        if dest_org_id:
            self.dest_org_id=dest_org_id
        else:
            raise Exception("dest_org_id is a required argument when creating a GeneMap")

        if dest_annot_type:
            self.dest_annot_type=dest_annot_type
        else:
            raise Exception("dest_annot_type is a required argument when creating a GeneMap")

        

    def map(self, item):
        # Detect which type of thing is being passed
        if isinstance(item, str):
            return 'Atoh1a'
        
        else:
            raise Exception("Error: unsupported object type ({0}) passed to map function".format(item_type))


"""
def remap_df_genes(orig_df: pd.DataFrame, orthomap_file: str):
    # Remap the passed-in Dataframe to have gene indexes from the orthologous mapping file.
    # Read HDF5 file using Pandas read_hdf
    orthomap_df = pd.read_hdf(orthomap_file)
    # Index -> gs1 / id2 / gs2
    orthomap_dict = orthomap_df.to_dict()["id2"]
    # NOTE: Not all genes can be mapped. Unmappable genes do not change in the original dataframe.
    return orig_df.rename(index=orthomap_dict)
"""
