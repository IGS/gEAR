"""
loaderutils.py - Helper functions for loading gene annotation data into the gEAR database.

Provides insert helpers for genes, gene symbols and gene URLs, plus lookup caches
keyed by Ensembl ID, gene ID or gene symbol.
"""

import math

def add_gene(curs, ensembl_id, organism_id, gene_symbol, biotype):
    """
    Insert a gene and its primary gene symbol.

    Returns:
        int: The ID of the newly inserted gene row.
    """
    add_g = ("INSERT INTO gene (ensembl_id, organism_id, biotype) "
             "VALUES (%s, %s, %s)")

    curs.execute( add_g, (ensembl_id, organism_id, biotype) )
    gene_id = curs.lastrowid

    add_gene_symbol(curs, gene_id, gene_symbol, 1)

    return gene_id

def add_gene_symbol(curs, gene_id, label, is_primary):
    """
    Insert a gene symbol for a gene and return the new row ID.
    """
    add_gs = "INSERT INTO gene_symbol (gene_id, label, is_primary) VALUES (%s, %s, %s)"
    curs.execute( add_gs, (gene_id, label, is_primary) )
    return curs.lastrowid

def cache_gene_aliases(curs, lower=False):
    """
    Build a lookup of primary and secondary gene symbols for each gene.

    Args:
        curs: Database cursor.
        lower (bool): If True, lowercase all gene symbols.

    Returns:
        dict: Gene ID -> {'p': primary symbol, 's': list of other symbols}.
    """
    query = "SELECT gene_id, label, is_primary FROM gene_symbol ORDER BY is_primary DESC"
    curs.execute(query)

    genes = dict()

    for (gene_id, gene_symbol, is_primary) in curs:
        if lower == True:
            gene_symbol = gene_symbol.lower()

        if gene_id not in genes:
            genes[gene_id] = {'p': gene_symbol, 's': []}
        else:
            genes[gene_id]['s'].append(gene_symbol)

    return genes

def cache_genes_by_ensembl_id(curs):
    """
    Return a dict mapping Ensembl ID to gene ID.
    """
    query = ("SELECT id, ensembl_id FROM gene")
    curs.execute(query)

    genes = dict()

    for (id, ensembl_id) in curs:
        genes[ensembl_id] = id

    return genes

def cache_gene_symbols_by_gene_id(curs, lower=False):
    """
    Return a dict mapping gene ID to primary gene symbol, optionally lowercased.
    """
    query = "SELECT gene_id, label FROM gene_symbol WHERE is_primary = 1"
    curs.execute(query)
    genes = dict()

    for (id, gene_symbol) in curs:
        if lower == True:
            gene_symbol = gene_symbol.lower()

        genes[id] = gene_symbol

    return genes

def cache_genes_by_primary_gene_symbol(curs, lower=False):
    """
    Return a dict mapping primary gene symbol to gene ID, optionally lowercased.
    """
    query = "SELECT gene_id, label FROM gene_symbol WHERE is_primary = 1"
    curs.execute(query)

    genes = dict()

    for (id, gene_symbol) in curs:
        if gene_symbol not in genes:
            if lower == True:
                gene_symbol = gene_symbol.lower()
                
            genes[gene_symbol] = id

    return genes

def cache_genes_by_gene_symbol(curs, lower=False):
    """
    Return a dict mapping every gene symbol (primary or alias) to gene ID.

    The first gene ID seen for a symbol is kept. If lower is True, symbols are lowercased.
    """
    query = "SELECT gene_id, label FROM gene_symbol"
    curs.execute(query)

    genes = dict()

    for (id, gene_symbol) in curs:
        if gene_symbol not in genes:
            if lower == True:
                gene_symbol = gene_symbol.lower()
                
            genes[gene_symbol] = id

    return genes

def fold_change_from_mean(num, mn, fold_changes):
    """
    Compute the log2 fold change of num relative to mean mn.

    The result is appended to fold_changes and also returned.
    """
    dev = num / mn
    fold_changes.append(math.log(dev, 2))
    return math.log(dev, 2)

def add_gene_url(curs, gene_id, label, url):
    """
    Insert a labeled URL for a gene and return the new row ID.
    """
    print("DEBUG: Inserting: {0} - {1}".format(label, url))
    qry = ("INSERT INTO gene_urls (gene_id, label, url) VALUES (%s, %s, %s)")
    curs.execute(qry, (gene_id, label, url))
    return curs.lastrowid


