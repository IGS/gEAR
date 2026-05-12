#!/opt/bin/python3

"""
This script takes as input a directory which contains these three files:

  expression.tab
  genes.tab
  observations.tab

And does the same process as the web-based uploader to convert it to an H5AD file.

"""

import argparse
import os
import re
import sys

import pandas as pd
import scanpy as sc
from scipy import sparse


def main():
    parser = argparse.ArgumentParser( description='3-tab -> H5AD')
    parser.add_argument('-i', '--input_directory', type=str, required=True, help='Path to an input directory containing the 3-tab files')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Path to the output file to be created')
    parser.add_argument('-r', '--row_chunk_size', type=int, required=False, default=500, help='Rows of input to be read at a time (helps control memory)')
    
    args = parser.parse_args()

    expression_matrix_path = None
    obs = None
    var = None

    # Get the base directory of the output file and make sure it's writable
    output_base_dir = os.path.dirname(args.output_file)
    if not os.path.exists(output_base_dir):
        print("Output directory does not exist: {0}".format(output_base_dir), file=sys.stderr, flush=True)
        sys.exit(1)
        
    if not os.access(output_base_dir, os.W_OK):
        print("Output directory is not writable: {0}".format(output_base_dir), file=sys.stderr, flush=True)
        sys.exit(1)        

    for infile in os.listdir(args.input_directory):
        # skip any files beginning with a dot
        if infile.startswith('.'):
            continue
        
        filepath = "{0}/{1}".format(args.input_directory, infile)
        
        # Read each file as pandas dataframes
        if infile == 'expression.tab' or os.path.basename(filepath)== 'expression.tab' or 'DataMTX.tab' in infile:
            expression_matrix_path = filepath
        elif infile == 'observations.tab' or os.path.basename(filepath)== 'observations.tab' or 'COLmeta.tab' in infile:
            print("Reading observations file: {0}".format(filepath), file=sys.stderr, flush=True)
            obs = pd.read_table(filepath, sep='\t', index_col=0, header=0, low_memory=False)
        elif infile == 'genes.tab' or os.path.basename(filepath)== 'genes.tab' or 'ROWmeta.tab' in infile:
            print("Reading genes file: {0}".format(filepath), file=sys.stderr, flush=True)
            var = pd.read_table(filepath, sep='\t', index_col=0, header=0)

    for str_type in ['cell_type', 'condition', 'time_point', 'time_unit']:
        if str_type in obs.columns:
            obs[str_type] = pd.Categorical(obs[str_type])

    for num_type in ['replicate', 'time_point_order']:
        if num_type in obs.columns:
            obs[num_type] = pd.to_numeric(obs[num_type])

    # Read in expressions as AnnData object in a memory-efficient manner
    print("Creating AnnData object with obs and var", file=sys.stderr, flush=True)
    adata = sc.AnnData(obs=var, var=obs)
    print("Reading expression matrix file: {0}".format(expression_matrix_path), file=sys.stderr, flush=True)    
    reader = pd.read_csv(expression_matrix_path, sep='\t', index_col=0, chunksize=args.row_chunk_size)
    
    ## Try to process the file the quickest way first, assuming things are peachy. Then, if not,
    #  do some checks and conversions (slower) as a backup.
    try:
        adata.X = sparse.vstack([sparse.csr_matrix(chunk.values) for chunk in reader])

    except Exception as e:
        print(f"\nOriginal vstack failed: {e}")
        print("Retrying with per-chunk cleanup...")

        expression_matrix = []
        chunk_shapes = []

        # Re-open reader here
        reader = pd.read_csv(expression_matrix_path, sep='\t', index_col=0, chunksize=args.row_chunk_size)

        for chunk_index, chunk in enumerate(reader, start=1):
            try:
                 # Clean each cell: strip string values
                chunk = chunk.replace(r'^\s+|\s+$', '', regex=True)
                chunk = chunk.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))

                # Convert to numeric (non-numeric → NaN → fill with 0)
                chunk_numeric = chunk.apply(pd.to_numeric, errors='coerce').fillna(0)

                matrix = sparse.csr_matrix(chunk_numeric.values)
                expression_matrix.append(matrix)
                chunk_shapes.append(matrix.shape)

            except Exception as inner_e:
                print(f"\nError in chunk {chunk_index}: {inner_e}")
                print("Chunk head:")
                print(chunk.head())
                raise

        # Try stacking the cleaned chunks
        try:
            adata.X = sparse.vstack(expression_matrix)
        except Exception as final_e:
            print(f"\nFinal vstack still failed: {final_e}")
            
            print("Collected chunk shapes:")
            for i, shape in enumerate(chunk_shapes):
                print(f"  Chunk {i+1}: {shape}")
                
            raise
        
    print("Finished reading expression matrix file", file=sys.stderr, flush=True)
    adata = adata.transpose()

    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'object':
            adata.obs[col] = adata.obs[col].fillna('').astype(str)
            
    adata.write(args.output_file)
            
if __name__ == '__main__':
    main()







