#!/opt/bin/python3

"""

"""

import cgi
import json
import os
import subprocess
import sys

def main():
    form = cgi.FieldStorage()
    #tmp_image_file = "/tmp/r_{0}.png".format(os.getpid())
    tmp_image_base = "../img/tmp_analysis/r_{0}_".format(os.getpid())
    
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
                                                                         method, tmp_image_base, cluster_count),
                            stdout=the_abyss,
                            shell=True)
    else:
        raise Exception("Other methods not yet handled");

    jdata = {
        'images': list()
    }

    for i in range(0, int(cluster_count)):
        web_path = './img/tmp_analysis/r_{0}_{1}.png'.format(os.getpid(), i + 1)
        jdata['images'].append({'path': web_path})

    print('Content-Type: application/json\n\n')
    print(json.dumps(jdata))

    
if __name__ == '__main__':
    main()
