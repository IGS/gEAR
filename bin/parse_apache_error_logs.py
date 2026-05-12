#!/usr/bin/env python3

"""

Look through the apache error logs and find/report errors for review.

Example log file names:
   ssl_umgear_error.log
   ssl_umgear_access.log.2.gz
   error.log
   error.log.1
   error.log.3.gz

Example log file lines:

My laptop:
[Wed Aug 21 23:22:11.960554 2019] [cgid:error] [pid 6597] [client ::1:52838] End of script output before headers: get_h5ad_dataset_list.cgi, referer: http://localhost/compare_datasets.html
[Wed Aug 21 23:04:50.670717 2019] [cgid:error] [pid 4028] [client ::1:43880] End of script output before headers: get_h5ad_data
set_list.cgi, referer: http://localhost/compare_datasets.html
Traceback (most recent call last):
  File "/var/www/cgi/get_h5ad_dataset_list.cgi", line 66, in <module>
    main()
  File "/var/www/cgi/get_h5ad_dataset_list.cgi", line 46, in main
    shared_with_user_collection.get_shared_with_user(has_h5ad=1, user=user, types=['microarray', 'bulk-rnaseq', 'singlecell-h5ad', 'single-cell-rnaseq', 'svg-expression', 'violin-standard'])
UnboundLocalError: local variable 'user' referenced before assignment

GCP:
[Sat Aug 24 14:15:50.408664 2019] [wsgi:error] [pid 16691] [remote 192.152.118.97:59224] saving figure to file /tmp/tsne_076a3a0c-c77e-4f36-aaa0-f37ff0ee80c8.png
[Sat Aug 24 14:16:03.447475 2019] [cgi:error] [pid 21624] [client 192.152.118.97:59223] AH01215: Traceback (most recent call last):: /var/www/cgi/get_user_layouts.cgi, referer: https://umgear.org/index.html
[Sat Aug 24 14:16:03.447575 2019] [cgi:error] [pid 21624] [client 192.152.118.97:59223] AH01215:   File "/var/www/cgi/get_user_layouts.cgi", line 55, in <module>: /var/www/cgi/get_user_layouts.cgi, referer: https://umgear.org/index.html
[Sat Aug 24 14:16:03.447584 2019] [cgi:error] [pid 21624] [client 192.152.118.97:59223] AH01215:     main(): /var/www/cgi/get_user_layouts.cgi, referer: https://umgear.org/index.html
[Sat Aug 24 14:16:03.447612 2019] [cgi:error] [pid 21624] [client 192.152.118.97:59223] AH01215:   File "/var/www/cgi/get_user_layouts.cgi", line 35, in main: /var/www/cgi/get_user_layouts.cgi, referer: https://umgear.org/index.html
[Sat Aug 24 14:16:03.447625 2019] [cgi:error] [pid 21624] [client 192.152.118.97:59223] AH01215:     params = [user.id]: /var/www/cgi/get_user_layouts.cgi, referer: https://umgear.org/index.html
[Sat Aug 24 14:16:03.447651 2019] [cgi:error] [pid 21624] [client 192.152.118.97:59223] AH01215: AttributeError: 'NoneType' object has no attribute 'id': /var/www/cgi/get_user_layouts.cgi, referer: https://umgear.org/index.html

"""

import argparse
import os


def main():
    parser = argparse.ArgumentParser( description='Apache error log parser')

    parser.add_argument('-i', '--input_dir', type=str, required=False, default="/var/log/apache2", help='Base apache log directory' )
    parser.add_argument('-m', '--mode', type=str, required=False, default='initcap', help='Mode of case in output: upper, lower, or initcap' )
    args = parser.parse_args()

    log_files = get_error_log_paths(args.input_dir)


    for file in log_files:
        print("Processing: {0}".format(file))



def get_error_log_paths(base):
    """
    Returns an array of files sorted by modified date, filtered for any which
    have 'error' in the file name.
    """
    paths = list()

    for filename in os.listdir(base):
        if 'error' in filename:
            paths.append("{0}/{1}".format(base, filename))

    paths.sort(key=lambda x: os.path.getmtime(x))

    return paths

if __name__ == '__main__':
    main()







