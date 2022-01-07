# Copying gene symbols between servers

## Server 1 (old server)

`mysql -u gear -h localhost -p -e "select g.ensembl_id, gs.label from gene g, gene_symbol gs where g.id = gs.gene_id" gear_portal > /tmp/gene_mapping.txt`

scp gene_mapping.txt to Server 2

## Server 2 (new server)

Highest valid `gene` id before updating - 290352

run `/opt/Python-3-7-3/bin/python3 <gear\_root>bin/load_old_gcid_gene_symbols.py /tmp/gene_mapping.txt