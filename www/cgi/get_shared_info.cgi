#!/opt/bin/python3

"""
Given both dataset and layout share IDs checks to see if there's a owner-curated
display with selected gene information. If so, returns the gene information.

Also returns some basic dataset information for display for a dataset share.

If only a layout share is given, the first dataset in the layout is retrieved 
and information is returned for that dataset.

Returns a JSON object with the following

{
    gene_symbol: 'symbol',
    dataset_label: 'label',
    layout_label: 'label',
    owner_name: 'Sally Nguyen'
}
"""

import cgi, json
import os, sys
import configparser
import json

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    config = configparser.ConfigParser()
    config.read('../../gear.ini')

    form = cgi.FieldStorage()
    dataset_share_id = form.getvalue('dataset_share_id')
    layout_share_id = form.getvalue('layout_share_id')
    dataset = None

    result = { 'gene_symbol': None, 'dataset_label': None, 
              'layout_label': None, 'owner_name': None }
    
    if dataset_share_id and dataset_share_id != 'null':
        dataset = geardb.get_dataset_by_share_id(share_id=dataset_share_id)

        if dataset:
            result['dataset_label'] = dataset.title
            owner = geardb.get_user_by_id(dataset.owner_id)

    elif layout_share_id and layout_share_id != 'null':
        layout = geardb.get_layout_by_share_id(layout_share_id)
        if layout:
            result['layout_label'] = layout.label
            owner = geardb.get_user_by_id(layout.user_id)
            layout.get_members()

            # just get the first dataset in the layout and use that one for search
            for layout_display in layout.members:
                dataset = geardb.get_dataset_by_id(id=layout_display.dataset_id)
                break

    if owner:
        result['owner_name'] = owner.user_name

    # Get the gene symbol for the dataset owner's default curation, if one exists.
    if dataset and owner:
        display = dataset.get_owner_display(is_multigene = False)
        if display:
            try:
                display_json = json.loads(display.plotly_config)
                if 'gene_symbol' in display_json:
                    result['gene_symbol'] = display_json['gene_symbol']

            except json.JSONDecodeError:
                display_json = {}
               
    print(json.dumps(result))

if __name__ == '__main__':
    main()
