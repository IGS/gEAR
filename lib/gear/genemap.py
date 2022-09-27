import os
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

      - A GeneCart object
      - A Dataset object
      - A single gene symbol string

    In each case, the source organism and dest organism IDs and annotation sources
    must also be defined.

    Example usage:

       gm = GeneMap(source_org_id=1, dest_org_id=1, source_annot_type='ensembl', dest_annot_type='ensembl')

       gc = GeneCart(...)
       filtered_gc = gm.map(gc)

    or:

       ds = Dataset(...)
       filtered_ds = gm.map(ds)

    or:

       mapped_gene = gm.map("Pou4f3")

    Any unmapped genes are simply excluded from the returned data structure, or "None" is
    returned in the case of an unmapped string.

    """
    def __init__(self, source_org_id=None, source_annot_type=None,
                 dest_org_id=None, dest_annot_type=None):
        self.source_org_id = source_org_id
        self.source_annot_type = source_annot_type
        self.dest_org_id=dest_org_id
        self.dest_annot_type=dest_annot_type

    def map(self, item):
        pass
