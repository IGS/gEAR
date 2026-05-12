# Copying gene symbols between servers

## Server 1 (old server)

`mysql -u gear -h localhost -p -e "select g.ensembl_id, gs.label from gene g, gene_symbol gs where g.id = gs.gene_id" gear_portal > /tmp/gene_mapping.txt`

scp gene_mapping.txt to Server 2

## Server 2 (new server)

Highest valid `gene` id before updating - 290352

run `/opt/Python-3-7-3/bin/python3 <gear\_root>bin/load_old_gcid_gene_symbols.py /tmp/gene_mapping.txt

## Misc GCID notes

In load_old_gcid_gene_symbols.py, I only added genes whose ensembl ID matched one already in the database.  This does create a bit of duplication in ensembl id entries, so I have kept track of periodic "last indexes" in case we want to reload.

Last gene in index before adding aspergillus genes - 326219

Last gene in index before adding rhizopus genes - 335835

* Ensembl IDs were not already in DB for Rhizopus so I have added them with a bunch of NULLs for fields
* I gave the Ensembl Release a value of -1 and so these never get used if the same gene has an earlier entry with a release value of 0 or higher.

Current last index - 340564