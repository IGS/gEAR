#!/usr/bin/env python3

"""
Changes the directory where schematic images are stored in schematic_filename
of dataset table.

Moves schematic file to new location.

"""

import os
import sys
import shutil

lib_path = os.path.abspath(os.path.join('..', 'lib'))
sys.path.append(lib_path)

import geardb


def main():

    schematic_dir = 'img/dataset_schematics'
    update_sql = "UPDATE dataset SET schematic_image = '{0}' WHERE id = '{1}';"
    update_list = list()
    count = 0

    cnx = geardb.Connection()
    cursor = geardb.Connection().get_cursor()

    qry = "SELECT id, schematic_image FROM dataset WHERE schematic_image IS NOT NULL"
    cursor.execute(qry)
    rows = cursor.fetchall()
    if rows is None:
        print("No datasets found with schematic images. Exiting.")
        return

    for (id, schematic_image) in rows:
        #skip gosling datasets and any ' ' values (not sure how that one happened)
        if 'gosling' not in schematic_image and len(schematic_image) > 2:
            print("Updating dataset: ", id)
            filename = schematic_image.rsplit('/', 1)[1]

            new_filepath = schematic_dir + "/" + filename

            #Add update sql statment
            update_list.append(update_sql.format(new_filepath, id))

            # Move schematic to new location
            old_path = '../www/' + schematic_image
            new_path = '../www/' + new_filepath

            #Move the schematic file if it's still in the old location
            if os.path.isfile(old_path):
                print("\tMoving file to: ", new_path)
                shutil.move(old_path, new_path)


    #Update the dataset table
    for query in update_list:
        cursor.execute(query)
        cnx.commit()

        count += 1

    cursor.close()
    print( "Finished.\n{0} datasets were updated.".format(str(count)) )



if __name__ == '__main__':
    main()
