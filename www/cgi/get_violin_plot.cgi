#!/opt/bin/python3

"""
Produces a violin plto using the seaborn library using an H5AD file as its source.

7f194001-38cd-77e4-7738-6575a9418350  gene:3871

sudo apt-get install gfortran libopenblas-dev liblapack-dev libpng-dev libfreetype6-dev
pip3 install seaborn

"""

import base64
import cgi
import json
import os
import sys

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDERR, killing this CGI.  So here we redirect
#  STDOUT until we need it.
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

import numpy as np
import scanpy as sc
sc.settings.verbosity = 0

def main():
    print('Content-Type: application/json\r\n\r\n')

    result = { 'success': 0, 'file': '' }

    form = cgi.FieldStorage()
    dataset_id = form.getvalue('dataset_id')
    gene_id = form.getvalue('gene_id')

    gene = geardb.get_gene_by_id(gene_id)
    ds = geardb.Dataset(id=dataset_id, has_h5ad=1)
    h5_path = ds.get_file_path()

    # Let's not fail if the file isn't there
    if not os.path.exists(h5_path):
        result['success'] = -1
        result['error'] = "no h5 file found for dataset {0}".format(dataset_id)
        print(json.dumps(result))
        sys.exit()

    headings = ['Expression', 'Location']

    adata = sc.read_h5ad(h5_path)
    # TODO: this is just for testing.  We need to do this on the datafiles as a pre-process
    adata.obs['mygroups'] = [name.split('--')[0] for name in adata.obs_names]

    image_dir = '/tmp/violins'
    pid = os.getpid()

    if not os.path.exists(image_dir):
        os.mkdir(image_dir)
    
    #sc.settings.figdir = image_dir # scanpy issue #73

    # seaborn.violinplot options: https://seaborn.pydata.org/generated/seaborn.violinplot.html
    try:
        ax = sc.pl.violin(adata, gene.gene_symbol, group_by='mygroups', show=False,
                          linewidth=0.5, scale="width", bw=.3)
    except IndexError:
        result['success'] = -1
        sys.stdout = original_stdout
        print('Content-Type: application/json\r\n\r\n')
        print(json.dumps(result))
        sys.exit(0)
        
    temp_image_file = "{0}/violin{1}.png".format(image_dir, pid)

    labels = []
    for l in ax.get_xticklabels():
        labels.append(l.get_text().replace('_', ' '))

    ax.set_xticklabels(labels, rotation=60)
    ax.set(ylabel='Expression')
    ax.set(xlabel='Location')
    ax.tick_params(axis='both', labelsize=12)

    fig = ax.get_figure()
    fig.tight_layout() # commenting this makes the graph fill better but truncates the bottom labels
    fig.savefig(temp_image_file, dpi=200)

    # convert image file to string for json (source: https://stackoverflow.com/a/35674047/2900840)
    with open(temp_image_file, 'rb') as f:
        result['file'] += base64.b64encode(f.read()).decode('ascii')

    result['success'] = 1
    os.remove(temp_image_file)

    sys.stdout = original_stdout
    print('Content-Type: application/json\r\n\r\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
