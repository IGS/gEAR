#!/opt/bin/python3

"""
Build an miRNA data list which consists of:
- miRNA Stem-loop and their mature forms (from miRBase miRNA.txt)
- miRNA aliases (from miRBase aliases.txt)
- miRNA families (from miRBase miFam.txt)

The idea is to construct this list; then use it for parsing EMBL or Genbank (not sure which yet) annotation files

Do coordinates match?
---------------------
Database	ID				Chr		Start   	Stop	Strand
Ensembl ENSGALG00000018303  24      3335097..3335176	+
miRBASE	MI0001258			24	    3330837..3330916	+

======> Ensembl vs miRBase start/stop vary, but Ensembl uses Gallus gallus 5.0 while miRBase is still 4.0.

Ensembl to miRBase relationship -- 1:1 or 1:Many?
-------------------------------------------------
ensembl_id          mirbase_id  ratio
ENSGALG00000025427	MI0007552	1:1
ENSDARG00000081501	MI0001958	1:1
ENSDARG00000081476	MI0001979	1:1
ENSDARG00000098981	MI0010848	1:1

======> Each Ensembl id is linked to 1 miRBase id.
"""


import argparse
import configparser
import cgi, json
import csv
import mysql.connector
from mysql.connector import errorcode
import re
import os
import sys

sys.path.append("{0}/../".format(os.path.dirname(sys.argv[0])))
import loaderutils

def main():
    parser = argparse.ArgumentParser( description='Builds list of miRNAs and their relationships')

    parser.add_argument('-a', '--alias_file', type=str, required=True, help='Path to  miRBase alias miRNA file (aliases.txt)' )
    parser.add_argument('-f', '--family_file', type=str, required=True, help='Path to miRBase family miRNA file (miFam.txt)' )
    parser.add_argument('-p', '--pub_file', type=str, required=True, help='Path to miRBase published miRNA file (miRNA.txt)' )
    parser.add_argument('-m', '--mgi_file', type=str, required=True, help='Path to MGI miRNA file (MGI_EntrezGene.rpt)' )
    args = parser.parse_args()

    # This will hold all miRNA info
    miRNAs = {}

    print("Parsing published data file")
    pub_file = open(args.pub_file)

    '''
        //
        ID   gga-mir-218-1     standard; RNA; GGA; 109 BP.
        XX
        AC   MI0001212;
        XX
        DE   Gallus gallus miR-218-1 stem-loop
        XX
        RN   [1]
        RX   PUBMED; 15592404.
        RA   International Chicken Genome Sequencing Consortium;
        RT   "Sequence and comparative analysis of the chicken genome provide unique
        RT   perspectives on vertebrate evolution";
        RL   Nature. 432:695-716(2004).
        XX
        RN   [2]
        RX   PUBMED; 23034410.
        RA   Meunier J, Lemoine F, Soumillon M, Liechti A, Weier M, Guschanski K, Hu H,
        RA   Khaitovich P, Kaessmann H;
        RT   "Birth and expression evolution of mammalian microRNA genes";
        RL   Genome Res. 23:34-45(2013).
        XX
        DR   ENTREZGENE; 777816; MIR218-1.
        XX
        FH   Key             Location/Qualifiers
        FH
        FT   miRNA           24..44
        FT                   /accession="MIMAT0001144"
        FT                   /product="gga-miR-218-5p"
        FT                   /evidence=experimental
        FT                   /experiment="Illumina [2]"
        FT                   /similarity="MI0000294"
        FT   miRNA           66..87
        FT                   /accession="MIMAT0026518"
        FT                   /product="gga-miR-218-3p"
        FT                   /evidence=experimental
        FT                   /experiment="Illumina [2]"
        XX
        SQ   Sequence 109 BP; 28 A; 19 C; 27 G; 0 T; 35 other;
             ugauaaugua gcgagauuuu cuguugugcu ugaucuaacc augugguugu gagguaugag        60
             uaaaacaugg uucugucaag caccauggaa cgucacgcag cuuucuaca                   109

    '''

    for line in pub_file:
        line.rstrip()
        cols = line.split()

        if line.startswith('//'):

            if stem_loop_label.startswith( ('hsa', 'gga', 'mmu', 'dre') ):
            # if stem_loop_label.startswith('gga'):
                # add the last mature miRNA in the stem-loop entry
                if mature_id not in qualifiers:
                    qualifiers[mature_id] = { 'primary_symbol_label': mature_primary_symbol,
                                              'similarity_id': similarity_id }
                if stem_loop_id not in miRNAs:
                    miRNAs[stem_loop_id] = { 'stem_loop_label': stem_loop_label,
                                             'mature_miRNAs': qualifiers }

            stem_loop_label = ''
            stem_loop_id = ''

        if line.startswith('ID'):
            stem_loop_label = cols[1]
            mature_miRNA_counter = 0

        if line.startswith('AC'):
            stem_loop_id = cols[1].strip(';')

        if line.startswith('FT'):
            if cols[1].startswith('miRNA'):
                if mature_miRNA_counter == 0:
                    qualifiers = {}
                else:
                    if mature_id not in qualifiers:
                        qualifiers[mature_id] = { 'primary_symbol_label': mature_primary_symbol,
                                                  'similarity_id': similarity_id }

                mature_id = ''
                mature_primary_symbol = ''
                similarity_id = ''

                mature_miRNA_counter += 1
            else:
                qualifier = cols[1].split('=')
                if qualifier[0] == '/accession':
                    mature_id = qualifier[-1].strip('"')

                if qualifier[0] == '/product':
                    mature_primary_symbol = qualifier[-1].strip('"')

                if qualifier[0] == '/similarity':
                    similarity_id = qualifier[-1].strip('"')

    # print(json.dumps(miRNAs, sort_keys=True, indent=4))
    '''
    The above code block yields:
        "MI0001212": {
            "mature_miRNAs": {
                "MIMAT0001144": {
                    "primary_symbol_label": "gga-miR-218-5p",
                    "similarity_id": "MI0000294"
                },
                "MIMAT0026518": {
                    "primary_symbol_label": "gga-miR-218-3p",
                    "similarity_id": ""
                }
            },
            "stem_loop_label": "gga-mir-218-1"
        },
    '''
    pub_file.close()

    print("Parsing alias file")
    # Parse alias file to get list of all stem-loop and mature miRNAs
    alias_file = open(args.alias_file)
    '''
    Looks like this:
        MI0000001	cel-let-7L;cel-let-7;
    '''

    for line in alias_file:
        line = line.rstrip()
        cols = line.split("\t")

        if cols[1].startswith( ('hsa', 'gga', 'mmu', 'dre') ):
        # if cols[1].startswith('gga'):
            miRNA_id = cols[0]
            aliases = cols[1]
            alias_list = aliases[:-1].split(";") #remove ';' ending and then split into a list

            # Identify primary symbol
            primary_sym = alias_list[-1]

            # Remove any duplicated aliases
            if len(alias_list) > 0:
                alias_list = list(set(alias_list))

            # Remove the primary_sym (if still present)
            if primary_sym in alias_list:
                alias_list.remove(primary_sym)

            for stem_loop_id in miRNAs:
                if miRNA_id.startswith('MIMAT'):
                    if miRNA_id in miRNAs[stem_loop_id]['mature_miRNAs']:
                    # if stem_loop_id == miRNA_id:
                        miRNAs[stem_loop_id]['mature_miRNAs'][miRNA_id]['aliases'] = alias_list

                else:
                    if stem_loop_id == miRNA_id:
                        miRNAs[stem_loop_id]['aliases'] = alias_list

    # print(json.dumps(miRNAs, sort_keys=True, indent=4))
    '''
    After aliases are added the same example is now structured like this:
        "MI0001212": {
            "aliases": [],
            "mature_miRNAs": {
                "MIMAT0001144": {
                    "aliases": [
                        "gga-miR-218"
                    ],
                    "primary_symbol_label": "gga-miR-218-5p",
                    "similarity_id": "MI0000294"
                },
                "MIMAT0026518": {
                    "aliases": [],
                    "primary_symbol_label": "gga-miR-218-3p",
                    "similarity_id": ""
                }
            },
            "stem_loop_label": "gga-mir-218-1"
        },
    '''
    alias_file.close()

    print("Parsing family file")
    #Parse family file to get stem_loop to family relationships
    family_file = open(args.family_file)
    '''
    Looks like this
        //
        AC   MIPF0000001
        ID   mir-17
        MI   MI0000071  hsa-mir-17
    '''

    stem_loop_id = None
    stem_loop_label = None
    family_id = None
    family_label = None

    for line in family_file:
        if line.startswith('//'):
            continue

        line = line.rstrip()
        cols = line.split()
        if cols[0] == 'AC':
            family_id = cols[1]

        if cols[0] == 'ID':
            family_label = cols[1]

        if cols[0] == 'MI':

            if cols[2].startswith( ('hsa', 'gga', 'mmu', 'dre') ):
            # if cols[2].startswith('gga'):
                stem_loop_id = cols[1]

                for miRNA in miRNAs:
                    if miRNA == stem_loop_id:
                        miRNAs[miRNA]['family_label'] = family_label
                        miRNAs[miRNA]['family_id'] = family_id

    # print(json.dumps(miRNAs, sort_keys=True, indent=4))
    '''
    After adding the family info, the same example from above is now:
        "MI0001212": {
            "aliases": [],
            "family_id": "MIPF0000026",
            "family_label": "mir-218",
            "mature_miRNAs": {
                "MIMAT0001144": {
                    "aliases": [
                        "gga-miR-218"
                    ],
                    "primary_symbol_label": "gga-miR-218-5p",
                    "similarity_id": "MI0000294"
                },
                "MIMAT0026518": {
                    "aliases": [],
                    "primary_symbol_label": "gga-miR-218-3p",
                    "similarity_id": ""
                }
            },
            "stem_loop_label": "gga-mir-218-1"
        },
    '''

    family_file.close()

    config = configparser.ConfigParser()
    config.read('../gear.ini')

    cnx = None
    try:
        cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                      host=config['database']['host'], database=config['database']['name'])
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password", file=sys.stderr)
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist", file=sys.stderr)
        else:
            print(err, file=sys.stderr)

    print("\nConnected to gEAR database.\n")

    cursor = cnx.cursor(buffered=True)

    # Get aliases from gEAR database
    cached_genes = loaderutils.cache_genes_by_gene_symbol(cursor, lower=False)
    cached_primary_genes = loaderutils.cache_genes_by_primary_gene_symbol(cursor, lower=False)
    print("gEAR genes have been cached.\n")

    '''
    Example of MGI file: (tab-delimited)

    MGI:2676793	Mirlet7a-1	O	microRNA let7a-1	24.84	13	Gene		387244	Mirnlet7a|mmu-let-7a-1|Mirnlet7a-1|Let-7a	miRNA gene	48538179	48538272	-	miRNA|ncRNA

    Source: http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt
    '''

    print("Parsing MGI_EntrezGene file...\n")
    mgi_file = open(args.mgi_file)
    reader = csv.reader(mgi_file, delimiter='\t')

    for cols in reader:
        if cols[0].startswith('MGI:'):
            mgi_gene_sym = cols[1]
            mgi_gene_name = cols[3] #column: Marker Name

            #Only interested in miR entries
            if mgi_gene_sym.startswith( ('mir', 'Mir') ) and mgi_gene_name.startswith('microRNA'):
                    # Skip MGI miR not found in gEAR
                    if mgi_gene_sym not in cached_genes and mgi_gene_sym.lower() not in cached_genes:
                        continue
                    else:
                        mgi_aliases = cols[9] #column: 'Synonyms |-delimited'
                        mgi_aliases_list = []
                        # Multiple aliases are '|' separated
                        if '|' in mgi_aliases:
                            mgi_aliases_list = mgi_aliases.split('|')
                        else:
                            mgi_aliases_list.append(mgi_aliases)

                        #Use the aliases_list to search the miRNAs object for a matching stem_loop_label
                        for alias in mgi_aliases_list:
                            for miRNA in miRNAs:
                                # If MGI alias matches stem_loop_label,
                                # change label to MGI symbol & then add MGI aliases to miRNAs object
                                if alias.lower() == miRNAs[miRNA]['stem_loop_label']:
                                    miRNAs[miRNA]['stem_loop_label'] = mgi_gene_sym

                                    if alias not in miRNAs[miRNA]['aliases']:
                                        miRNAs[miRNA]['aliases'].append(alias.lower())

    # Set up SQL statements
    update_sym_as_primary = '''
                            UPDATE gene_symbol
                            SET is_primary = 1
                            WHERE gene_id = %s
                                AND label = %s
                            '''

    add_miRNA_to_gene = '''
                        INSERT INTO gene
                        (ensembl_id, genbank_acc, organism_id, molecule, start, stop, gene_symbol, product, biotype, shield_facs_chart_url)
                        VALUES (NULL, NULL, %s, NULL, NULL, NULL, NULL, %s, 'miRNA', NULL)
                        '''

    add_miRNA_to_gene_symbol = '''
                        INSERT INTO gene_symbol
                        (gene_id, label, is_primary)
                        VALUES (%s, %s, %s)
                        '''

    miRNA_counter = 0
    miRNAs_added = {} # will hold labels that have been just added. Prevents duplicate entries
    for stem_loop_id in miRNAs:
        # NOTE: Should not have to check if stem-loops need to be added to gEAR database
        # They are added when annotation is updated.

        # Get the stem-loop's gene_id
        stem_loop_label = miRNAs[stem_loop_id]['stem_loop_label']

        if stem_loop_label not in cached_genes and stem_loop_label.lower() not in cached_genes:
            continue # skips those not in gEAR because they lack a chromosome location
        else:
            try:
                gene_id = cached_genes[stem_loop_label]
            except KeyError:
                gene_id = cached_genes[stem_loop_label.lower()]

            if stem_loop_label not in cached_primary_genes and stem_loop_label.lower() not in cached_primary_genes:
                #A primary sym is mistakenly marked as secondary...fix it
                cursor.execute(update_sym_as_primary, (gene_id, stem_loop_label))
                cnx.commit()
                print("{0} now set as primary gene symbol.".format(stem_loop_label))

            miRNAs[stem_loop_id]['gene_id'] = gene_id

        #Determine organism id
        if stem_loop_label.startswith( ('mmu', 'Mir', 'mir') ):
            organism_id = 1 # Mouse

        if stem_loop_label.startswith('hsa'):
            organism_id = 2 # Human

        if stem_loop_label.startswith('dre'):
            organism_id = 3 # Zebrafish

        if stem_loop_label.startswith('gga'):
            organism_id = 5 # Chicken

        # Check if stem-loop aliases need to be added to gene table
        stem_loop_aliases = miRNAs[stem_loop_id]['aliases']
        if len(stem_loop_aliases) > 0:
            miRNAs[stem_loop_id]['alias_ids'] = {}
            for stem_loop_alias in stem_loop_aliases:
                if stem_loop_alias not in cached_genes:
                    #Need to add alias symbol

                    #Insert stem-loop alias into 'gene_symbol' table
                    cursor.execute(add_miRNA_to_gene_symbol, (gene_id, stem_loop_alias, 0))
                    cnx.commit()
                    miRNA_counter += 1

                    #Prevent duplicates
                    # miRNAs_added[stem_loop_alias] = gene_id
                else:
                    # Already in database, so get the gene_id
                    gene_id = cached_genes[stem_loop_alias]
                    # print('Found alias miRNA: {0}. gene_id = {1}'.format(alias, gene_id))

                # save the mature miRNA's gene_id
                miRNAs[stem_loop_id]['alias_ids'][stem_loop_alias] = gene_id

        # Check if mature miRNAs need to be added to gEAR database (gene table)
        for mature_miRNA in miRNAs[stem_loop_id]['mature_miRNAs']:
            mature_primary_sym = miRNAs[stem_loop_id]['mature_miRNAs'][mature_miRNA]['primary_symbol_label']
            mature_aliases = miRNAs[stem_loop_id]['mature_miRNAs'][mature_miRNA]['aliases']

            if mature_primary_sym in miRNAs_added:
                # Already in database, so get the gene_id
                gene_id = miRNAs_added[mature_primary_sym]

            else:
                if mature_primary_sym not in cached_genes:
                    # Mature primary sym Needs to be added

                    #Insert mature miRNA into 'gene' table
                    cursor.execute(add_miRNA_to_gene, (organism_id, mature_primary_sym))
                    gene_id = cursor.lastrowid
                    cnx.commit()
                    # print("Added mature miRNA: {0}. gene_id = {1}".format(mature_primary_sym, gene_id))

                    #Insert primary_symbol into 'gene_symbol' table
                    cursor.execute(add_miRNA_to_gene_symbol, (gene_id, mature_primary_sym, 1))
                    cnx.commit()
                    miRNA_counter += 1

                    #Prevent duplicates
                    miRNAs_added[mature_primary_sym] = gene_id

                else:
                    # Redundant. Remove it?
                    # Already in database, so get the gene_id
                    gene_id = cached_genes[mature_primary_sym]
                    # print('Found mature miRNA: {0}. gene_id = {1}'.format(mature_primary_sym, gene_id))

            # save the mature miRNA's gene_id
            miRNAs[stem_loop_id]['mature_miRNAs'][mature_miRNA]['gene_id'] = gene_id

            # Now check if mature aliases need added to 'gene' table
            if len(mature_aliases) > 0:
                miRNAs[stem_loop_id]['mature_miRNAs'][mature_miRNA]['alias_ids'] = {}
                for mature_alias in mature_aliases:
                    if mature_alias in miRNAs_added:
                        # Already in database, so get the gene_id
                        gene_id = miRNAs_added[mature_primary_sym]
                    else:
                        if mature_alias not in cached_genes:
                            #Need to add alias symbol

                            #Insert stem-loop alias into 'gene_symbol' table
                            cursor.execute(add_miRNA_to_gene_symbol, (gene_id, mature_alias, 0))
                            cnx.commit()
                            miRNA_counter += 1

                            #Prevent duplicates
                            miRNAs_added[mature_alias] = gene_id
                        else:
                            # Already in database, so get the gene_id
                            gene_id = cached_genes[mature_alias]
                            # print('Found alias miRNA: {0}. gene_id = {1}'.format(alias, gene_id))

                    # save the mature miRNA's gene_id
                    miRNAs[stem_loop_id]['mature_miRNAs'][mature_miRNA]['alias_ids'][mature_alias] = gene_id

    # print('{0} mature miRNAs need to be added to gEAR "gene" table'.format(miRNA_counter))
    # print(json.dumps(miRNAs, sort_keys=True, indent=4))
    '''
    Continuing the example. miRNAs now includes 'gene_id' and looks like this:
        "MI0001212": {
            "aliases": [],
            "family_id": "MIPF0000026",
            "family_label": "mir-218",
            "gene_id": 210901,
            "mature_miRNAs": {
                "MIMAT0001144": {
                    "alias_ids": {
                        "gga-miR-218": 221150
                    },
                    "aliases": [
                        "gga-miR-218"
                    ],
                    "gene_id": 221149,
                    "primary_symbol_label": "gga-miR-218-5p",
                    "similarity_id": "MI0000294"
                },
                "MIMAT0026518": {
                    "aliases": [],
                    "gene_id": 221148,
                    "primary_symbol_label": "gga-miR-218-3p",
                    "similarity_id": ""
                }
            },
            "stem_loop_label": "gga-mir-218-1"
        },
    '''
    print("\n{0} miRNA gene symbols added to 'gene_symbol'\n".format(miRNA_counter))

    print("Checking if new miRNA relationships need added to mirna_family...\n")
    add_mirna_to_mirna_family = """
            INSERT INTO mirna_family
            (stem_loop_id, mature_id, family_id, family_label)
            VALUES (%s, %s,%s, %s)
    """
    # 1) insert stem-loop, mature_miRNAs, and family id & label (if present)
    # 2) (if stem-loop alias_ids present) insert stem-loop alias, mature_miRNAs, and family id & label (if present)
    # 3) (if mature_miRNA alias_ids present) insert stem-loop, mature_miRNAs, and family id & label (if present)
    # 4) (if stem-loop alias_ids present AND mature_miRNA alias_ids present) insert stem-loop, mature_miRNAs, and family id & label (if present)

    relationships_added = 0
    for stem_loop_id in miRNAs:

        # TODO: Added this until can address stem-loops that lack a chromosome assignment
        #       No chr location = No ensembl entry in annotation file ==> not added to gEAR
        #       Example: gga-mir-7480-1 MI0024159
        if 'gene_id' in miRNAs[stem_loop_id]:
            stem_loop_gene_id = miRNAs[stem_loop_id]['gene_id']

            # Is stem-loop in a family? (not all stem-loops belong to a family)
            if 'family_id' in miRNAs[stem_loop_id]:
                family_id = miRNAs[stem_loop_id]['family_id']
                family_label = miRNAs[stem_loop_id]['family_label']
            else:
                family_id = None
                family_label = None

            for mature_id in miRNAs[stem_loop_id]['mature_miRNAs']:
                mature_gene_id = miRNAs[stem_loop_id]['mature_miRNAs'][mature_id]['gene_id']

                # If not in database, add to mirna_family: stem-loop, mature, family id, family label
                if check_miRNA_family(cursor, stem_loop_gene_id, mature_gene_id, family_id, family_label) is False:
                    cursor.execute(add_mirna_to_mirna_family, (stem_loop_gene_id, mature_gene_id, family_id, family_label))
                    cnx.commit()
                    relationships_added += 1

                # If mature miRNA has any aliases with gene_ids, add them
                if 'alias_ids' in miRNAs[stem_loop_id]['mature_miRNAs'][mature_id]:
                    for mature_alias in miRNAs[stem_loop_id]['mature_miRNAs'][mature_id]['alias_ids']:
                        mature_alias_gene_id = miRNAs[stem_loop_id]['mature_miRNAs'][mature_id]['alias_ids'][mature_alias]

                        # If not in database, add to mirna_family: stem-loop, mature alias, family id, family label
                        if check_miRNA_family(cursor, stem_loop_gene_id, mature_alias_gene_id, family_id, family_label) is False:
                            cursor.execute(add_mirna_to_mirna_family, (stem_loop_gene_id, mature_alias_gene_id, family_id, family_label))
                            cnx.commit()
                            relationships_added += 1

            # If stem-loop has any aliases with gene_ids, add them
            if 'alias_ids' in miRNAs[stem_loop_id]:
                for alias in miRNAs[stem_loop_id]['alias_ids']:
                    stem_loop_alias_gene_id = miRNAs[stem_loop_id]['alias_ids'][alias]

                    for mature_id in miRNAs[stem_loop_id]['mature_miRNAs']:
                        mature_gene_id = miRNAs[stem_loop_id]['mature_miRNAs'][mature_id]['gene_id']

                        # If not in database, add to mirna_family: stem-loop alias, mature, family id, family label
                        if check_miRNA_family(cursor, stem_loop_alias_gene_id, mature_gene_id, family_id, family_label) is False:
                            cursor.execute(add_mirna_to_mirna_family, (stem_loop_alias_gene_id, mature_gene_id, family_id, family_label))
                            cnx.commit()
                            relationships_added += 1

                        # If mature miRNA has any aliases with gene_ids, add them
                        if 'alias_ids' in miRNAs[stem_loop_id]['mature_miRNAs'][mature_id]:
                            for mature_alias in miRNAs[stem_loop_id]['mature_miRNAs'][mature_id]['alias_ids']:
                                mature_alias_gene_id = miRNAs[stem_loop_id]['mature_miRNAs'][mature_id]['alias_ids'][mature_alias]

                                # If not in database, add to mirna_family: stem-loop alias, mature alias, family id, family label
                                if check_miRNA_family(cursor, stem_loop_alias_gene_id, mature_alias_gene_id, family_id, family_label) is False:
                                    cursor.execute(add_mirna_to_mirna_family, (stem_loop_alias_gene_id, mature_alias_gene_id, family_id, family_label))
                                    cnx.commit()
                                    relationships_added += 1


    print("{0} relationships added to 'mirna_family'\n".format(relationships_added))
    # print("Task Complete. miRNA relationships generated and inserted into 'mirna_family' table.")

    '''
    The relationships added to gEAR for the continuing example:

        mysql> select * from mirna_family m JOIN gene_symbol g ON g.gene_id=m.mature_id where m.stem_loop_id='210901';
        +-----+--------------+-----------+-------------+--------------+--------+---------+----------------+------------+
        | id  | stem_loop_id | mature_id | family_id   | family_label | id     | gene_id | label          | is_primary |
        +-----+--------------+-----------+-------------+--------------+--------+---------+----------------+------------+
        | 395 |       210901 |    221149 | MIPF0000026 | mir-218      | 261920 |  221149 | gga-miR-218-5p |          1 |
        | 396 |       210901 |    221150 | MIPF0000026 | mir-218      | 261921 |  221150 | gga-miR-218    |          0 |
        | 397 |       210901 |    221148 | MIPF0000026 | mir-218      | 261919 |  221148 | gga-miR-218-3p |          1 |
        +-----+--------------+-----------+-------------+--------------+--------+---------+----------------+------------+


    '''

    print("Checking if new gene_urls need to be added...\n")
    # Add gene_urls for miRBase
    # Get gene_urls
    gene_urls_sql = '''
                    SELECT gene_id, label
                    FROM gene_urls
                    '''
    cached_gene_urls = {}

    cursor.execute(gene_urls_sql,)
    for row in cursor:
        gene_id = row[0]
        label = row[1]

        if gene_id not in cached_gene_urls:
            cached_gene_urls[gene_id] = []

        if label not in cached_gene_urls[gene_id]:
            cached_gene_urls[gene_id].append(label)

    add_gene_url_sql = '''
                        INSERT INTO gene_urls
                        (gene_id, label, url)
                        VALUES (%s, %s, %s)
                        '''
    stem_loop_label = 'miRBase - Stem-Loop'
    mature_label = 'miRBase - Mature'
    stem_loop_url_base = 'http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc='
    mature_url_base = 'http://www.mirbase.org/cgi-bin/mature.pl?mature_acc='

    counter_new_urls = 0
    for stem_loop_id in miRNAs:
        if 'gene_id' in miRNAs[stem_loop_id]:
            stem_loop_gene_id = miRNAs[stem_loop_id]['gene_id']

            # Gene found in gene_urls, but no url for miRBase
            if stem_loop_gene_id in cached_gene_urls:
                if stem_loop_label not in cached_gene_urls[stem_loop_gene_id]:
                    stem_loop_url = stem_loop_url_base + stem_loop_id
                    cursor.execute(add_gene_url_sql, (stem_loop_gene_id, stem_loop_label, stem_loop_url))
                    cnx.commit()
                    counter_new_urls += 1

            # Gene not found in gene_urls, so no url for miRBase
            if stem_loop_gene_id not in cached_gene_urls:
                stem_loop_url = stem_loop_url_base + stem_loop_id
                cursor.execute(add_gene_url_sql, (stem_loop_gene_id, stem_loop_label, stem_loop_url))
                cnx.commit()
                counter_new_urls += 1

            for mature_id in miRNAs[stem_loop_id]['mature_miRNAs']:
                mature_gene_id = miRNAs[stem_loop_id]['mature_miRNAs'][mature_id]['gene_id']

                # Gene found in gene_urls, but no url for miRBase
                if mature_gene_id in cached_gene_urls:
                    if mature_label not in cached_gene_urls[mature_gene_id]:
                        mature_url = mature_url_base + mature_id
                        cursor.execute(add_gene_url_sql, (mature_gene_id, mature_label, mature_url))
                        cnx.commit()
                        counter_new_urls += 1

                # Gene not found in gene_urls, so no url for miRBase
                if mature_gene_id not in cached_gene_urls:
                    mature_url = mature_url_base + mature_id
                    cursor.execute(add_gene_url_sql, (mature_gene_id, mature_label, mature_url))
                    cnx.commit()
                    counter_new_urls += 1

    print("{0} new entries added to gene_urls.\n".format(counter_new_urls))

    print('\nFinished.\n')
    cursor.close()
    cnx.close()

def check_miRNA_family(cursor, stem_loop_id, mature_id, family_id, family_label):
    # Checks if the miRNA family--->stem-loop--->mature miR relationship already exists
    qry = """
        SELECT COUNT(id)
        FROM mirna_family
        WHERE stem_loop_id = {0}
          AND mature_id = {1}
    """.format(stem_loop_id, mature_id)

    cursor.execute(qry)
    for row in cursor:
        if row[0] > 0:
            # miRNA family entry already exists
            # print('{0} miR Fam entries found.'.format(row[0]))
            return True
        else:
            # No miRNA family entry found
            # print("Will add: stem-loop: {0}, mature: {1}".format(stem_loop_id, mature_id))
            return False

if __name__ == '__main__':
    main()
