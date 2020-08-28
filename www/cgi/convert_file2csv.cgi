#!/opt/bin/python3

"""
Accepts spreadsheet files and converts them to 'tab' delimited text files.
Not intended for .tab, .csv, .tsv, or .txt, files.

Uses Gnumeric's ssconvert to make the format conversion

command used to test:
ssconvert --export-type=Gnumeric_stf:stf_assistant -O "separator='  ' charset=ascii format=raw" ./uploads/files/example_bargraph_standard.ods ./uploads/files/example_bargraph_standard.tab

UPDATE: Now separates columns by ',' so files are truly converted into CSV files, not TAB as stated above. This allows for spaces in column heading names across formats
"""

import cgi, json
import subprocess
import sys
import os

sys.path.append('../..')
import loaderutils

def main():

    print('Content-Type: application/json\n\n')

    form = cgi.FieldStorage()
    input_file = form.getvalue('filename')
    filetype = form.getvalue('filetype')

    #ssconvert cannot handle spaced filenames. replace spaces with '_'
    os.rename('../uploads/files/{0}'.format(input_file), '../uploads/files/{0}'.format(input_file.replace(" ", "_")))

    #replace spaces with "_"
    input_file = input_file.replace(" ", "_")

    # get the file's output name
    output_file = input_file.replace(filetype, 'csv')

    # will be de-limiter for output file
    # separator = "'  '"
    separator = ","

    result = {'success': False, 'filename': ''}

    # export as csv, separate with 2 spaces, preserve format to keep decimal places
    # quoting-mode always - older version of ssconvert will drop values of 0, so add "" as placeholder
    convert_file_command = 'ssconvert --export-type=Gnumeric_stf:stf_assistant -O "separator={0} charset=ascii format=preserve quoting-mode=always" ../uploads/files/{1} ../uploads/files/{2}'.format(separator, input_file, output_file)

    #convert file to tab file
    subprocess.call(convert_file_command, shell=True)

    # check that conversion was successful
    result['success'] = True
    result['filename'] = output_file

    #remove original file
    os.remove('../uploads/files/{0}'.format(input_file))

    print(json.dumps(result))

if __name__ == '__main__':
    main()
