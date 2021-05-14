#!/opt/bin/python3

"""
We often want to generate a plot with R and make it the source for an image in gEAR,
but R does not have support for printing binary streams to the STDOUT buffer.  To
avoid writing tmp files everywhere, this script does the following:

1. Takes arguments to generate a plot with R
2. Save the result file in a temporary place
3. Print the image as a binary stream (PNG)
4. Delete the temp file

"""

import cgi
import os
import subprocess
import sys

def main():
    form = cgi.FieldStorage()
    tmp_image_file = "/tmp/r_{0}.png".format(os.getpid())
    gene = form.getvalue('gene')   # example: "E2F1"
    dataset_base = form.getvalue('dataset_id') # example: "BrainSpanDFC"
    dataset_path = "/usr/local/projects/gEAR/carlo/"

    # 'single_gene', 'pca', 'kmeans', 'hclust'
    method = form.getvalue('method');
    cluster_count = form.getvalue('cluster_count');

    if method == 'single_gene':
        with open(os.devnull, "w") as the_abyss:
            subprocess.call("./plotCode1gene.R {0} {1} {2} {3}".format(gene, dataset_base,
                                                                       dataset_path, tmp_image_file),
                            stdout=the_abyss,
                            shell=True)
    elif method == 'pca' or method == 'kmeans' or method == 'hclust':
        with open(os.devnull, "w") as the_abyss:
            subprocess.call("./plot_r_analysis.R {0} {1} {2} {3} {4}".format(dataset_base, dataset_path,
                                                                         method, tmp_image_file, cluster_count),
                            stdout=the_abyss,
                            shell=True)
    else:
        raise Exception("Other methods not yet handled");

    # write the PNG as a binary stream
    with open(tmp_image_file, 'rb') as f:
        print("Content-Type: image/png\n")
        sys.stdout.flush() # <---
        sys.stdout.buffer.write(f.read())

    # delete the tmp file
    os.remove(tmp_image_file)


if __name__ == '__main__':
    main()
