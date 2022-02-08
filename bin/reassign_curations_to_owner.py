#!/opt/bin/python3

"""
Not sure how generally useful this will be.

During a conference Bea did multi-gene curations for a LOT of datasets, and this 
reassigns ownership of those curations to to the dataset owner.

It also sets these curations as the default multi-gene view for each, clearing any
previous defaults.
"""
import argparse
import os
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import geardb

CURATOR_USER_ID = 499

def main():
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    dataset_ids = list()

    layout_names = [
        #'Adult',
        #'Aging',
        #'Central auditory system',
        #'Chick',
        #'Development',
        #'Ear (diverse variety)',
        #'Hair cell',
        #'Hair cell and supporting cell transcriptomics and epigenomics, Mouse (Segil 2021)',
        #'Hair cell reprogramming - transcriptomics and epigenetics (Segil)',
        #'Hearing (default)',
        #'Human disease',
        #'Lateral wall',
        #'Mouse models',
        #'Mouse models 2',
        #'Neuron'
    ]

    layout_qry = """
          SELECT id, user_id, label
            FROM layout 
           WHERE user_id = 499
    """
    qry_args = []

    cursor.execute(layout_qry, qry_args)

    for (layout_id, user_id, label) in cursor:
        if label in layout_names:        
            layout = geardb.Layout(id=layout_id)
            dataset_ids.extend(layout.dataset_ids())

    ## For each datasets, grab the displays
     # If a display is owned by curator change it to be owned by the dataset author
     # Also clear any preferences for that dataset for the author and set the new one
    dc = geardb.DatasetCollection()
    dc.get_by_dataset_ids(dataset_ids)

    for dataset in dc.datasets:
        displays = dataset.get_displays()

        for display in displays:
            if display.user_id == CURATOR_USER_ID and dataset.owner_id != CURATOR_USER_ID:
                ## change display to be owned by dataset owner
                display.user_id = dataset.owner_id
                #display.save()

                ## clear any existing preferences for this dataset for this user
                ## set this one as the new preference for multi-genes
                query = """
                        INSERT INTO dataset_preference (user_id, dataset_id, display_id, is_multigene)
                             VALUES (%s, %s, %s, %s)
                             ON DUPLICATE KEY UPDATE display_id = VALUES(display_id)
                """
                print("Executing update with ({0}, {1}, {2}, {3})".format(display.user_id, dataset.id, display.id, 1))
                #cursor.execute(query, (display.user_id, dataset.id, display.id, 1))
                
    cnx.commit()
    cursor.close()
    cnx.close()



    
    
if __name__ == '__main__':
    main()
    
