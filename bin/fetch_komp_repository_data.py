#!/usr/bin/env python3

"""
Link page with all genes:

https://www.komp.org/catalog.php?available=&mutation=&gname=all&project=&origin=&

From this, we find what we need by parsing (extra elements/attributes removed):

<table id='catalog'>
   <tbody>
      <tr>
         <td><a href='geneinfo.php?geneid=50'><b><i>0610007P14Rik</i></b></a></td>
	 <td>KOMP</td>
	 <td>CSD</td>
         <td><a href='geneinfo.php?geneid=50'>Products available.<br/>(ES Cells)</a></td>
      </tr>
   </tbody>
</table>

===========================================================================================
===========================================================================================

Synonyms are stored like:

<font style="font-size:90%"><i>Synonyms: </i>2510005N23Rik, 9930116O05Rik, D130086K05Rik, Rfxdc2</font>

"""

import argparse
from bs4 import BeautifulSoup
import pickle
import os
import re
import urllib.request

def main():
    parser = argparse.ArgumentParser( description='')

    SOURCE_URL = 'https://www.komp.org/catalog.php?available=&mutation=&gname=all&project=&origin=&'
    KOMP_BASE_URL = 'https://www.komp.org'

    ## output file to be written
    ## parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to an input file to be read' )
    args = parser.parse_args()

    #print("INFO: downloading source page ... ")
    #urllib.request.urlretrieve(SOURCE_URL, 'tmp.current_master.html')
    #print("done.")

    print("INFO: parsing source page ... ")
    soup = BeautifulSoup(open('tmp.current_master.html'))
    catalog = soup.find(id='catalog')

    # The catalog has duplicate entries of the same gene from different origins.
    #  This keeps track of them so I don't process the same ones twice.
    genes_found = dict()
    url_labels_to_skip = ['Home', 'MyKOMP', 'Log In', 'New Account', 'Catalog', 'Google', 'Make my mouse...!']

    dstore = dict()

    for row in catalog.tbody.find_all('tr'):
        gene_id = str(row.td.a.i.string)
        gene_url = row.td.a['href']

        if gene_id not in genes_found:
            genes_found[gene_id] = gene_url

    whitespace_pattern = re.compile(r'\s+')

    for gene_id in sorted(genes_found):
        print("Processing: {0} - {1}".format(gene_id, genes_found[gene_id]))
        urllib.request.urlretrieve("{0}/{1}".format(KOMP_BASE_URL, genes_found[gene_id]), 'tmp.gene.html')

        dstore[gene_id] = { 'syn': list(), 'links': list() }

        gene_soup = BeautifulSoup(open('tmp.gene.html'))

        font_tags = gene_soup.find_all('font')
        
        for font_tag in font_tags:
            #print("INFO: Font tag: {0}".format(font_tag.get_text() ) )

            m = re.search('Synonyms: (.+)', font_tag.get_text() )
            if m:
                #print("INFO: Synonyms: {0}".format(m.group(1)))
                dstore[gene_id]['syn'].extend( str(m.group(1)).split(',') )
                
        # All the links we care about here are nested within li elements, so grab them that way
        list_items = gene_soup.find_all('li')
        links = list()
        urls_found = dict()

        for li in list_items:
            links.extend( li.find_all('a') )

        for link in links:
            label = link.string

            if label is None or label in url_labels_to_skip:
                continue
            else:
                label = label.strip()

            href = re.sub(whitespace_pattern, '', link['href'])

            if 'target' in link:
                target = link['target']
            else:
                target = None

            if href in urls_found:
                continue
            else:
                urls_found[href] = 1
                
            #print("\t{0} - {1} - {2}".format(label, target, href))
            dstore[gene_id]['links'].append( {'label': str(label), 'href': str(href) } )
        
        # just limit to 1 for now
        #break

    pickle.dump( dstore, open("komp_repository_data.p", "wb") )
    
    
    


if __name__ == '__main__':
    main()







