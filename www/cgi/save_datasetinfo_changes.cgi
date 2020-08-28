#!/opt/bin/python3

"""
Changes the information of a dataset via .editable elements in the dataset_manager

This script first checks if the user owns the dataset, then proceeds with the info change
if they own it.

If the user does not own the dataset, an error is returned stating that.
Successful access_level change returns the dataset JSON that was changed.

Requires:
1) Session id - which contains user_id
2) Dataset id
3) One of the following:
    - Title
    - Description
    - Pubmed ID
    - GEO ID
    - Schematic File name

"""

import os
import cgi
import json
import os
import sys
import re
import shutil

lib_path = os.path.abspath(os.path.join('..', '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    print('Content-Type: application/json\n\n')

    cnx = geardb.Connection()
    cursor = cnx.get_cursor()
    form = cgi.FieldStorage()
    session_id = form.getvalue('session_id')
    dataset_id = form.getvalue('id')
    field2change = form.getvalue('field')
    new_value = form.getvalue('new_value')

    user = geardb.get_user_from_session_id(session_id)

    # Does user own the dataset...
    owns_dataset = check_dataset_ownership(cursor, user.id, dataset_id)

    if owns_dataset == True:
        result = { 'dataset':[] }

        if field2change == "title" or field2change == "ldesc" or field2change == "pubmed" or field2change == "geo":
            change_dataset_info(cursor, field2change, user.id, dataset_id, new_value)
            cnx.commit()

            # get the updated value...
            result['dataset'].extend(get_dataset(cursor, field2change, user.id, dataset_id))
            result['success'] = True

            cursor.close()
            cnx.close()

            print(json.dumps(result))

        elif field2change == 'tag':

            #Get list of tags already in gEAR
            qry_get_tags = """
                SELECT label, id
                FROM tag
            """
            cached_tags = {}
            cursor.execute(qry_get_tags)
            for row in cursor:
                cached_tags[row[0].lower()] = row[1]

            if new_value is None:
                #Remove all tags from dataset
                remove_tags_from_dataset_tag = """
                        DELETE FROM dataset_tag
                        WHERE dataset_id=%s;
                """
                cursor.execute(remove_tags_from_dataset_tag, (dataset_id,))
                cnx.commit()

                # get the updated value...
                result['dataset'].extend(get_dataset(cursor, field2change, user.id, dataset_id))
                result['success'] = True

                cursor.close()
                cnx.close()

                print(json.dumps(result))

            else:
                new_tags_list = []
                raw_tags = new_value.split(',')
                print(raw_tags, file=sys.stderr)

                # Ensures duplicates are removed
                for tag in raw_tags:
                    if tag not in new_tags_list:
                        new_tags_list.append(tag)

                qry_get_datasettags = """
                      SELECT t.label
                      FROM tag t
                      JOIN dataset_tag d ON d.tag_id=t.id
                      WHERE d.dataset_id=%s
                """
                #get old tags linked to dataset
                cursor.execute(qry_get_datasettags, (dataset_id,))

                #were old tags kept or removed?
                for row in cursor:
                    old_tag = row[0]

                    #old tag was removed
                    if old_tag not in new_tags_list:
                        #remove old_tag from dataset_tag
                        remove_tag_from_dataset_tag = """
                                DELETE d FROM dataset_tag d
                                LEFT JOIN tag t ON t.id=d.tag_id
                                WHERE t.label=%s;
                        """
                        cursor.execute(remove_tag_from_dataset_tag, (old_tag,))
                        cnx.commit()

                    #old tag was kept
                    else:
                        #remove old_tag new_tags_list. No need to add it again
                        old_tag_index = new_tags_list.index(old_tag)
                        del new_tags_list[old_tag_index]

                add_tag_sql = """
                    INSERT INTO tag (label)
                    VALUES (%s)
                """
                add_datasettag_sql = """
                    INSERT INTO dataset_tag (tag_id, dataset_id)
                    VALUES (%s, %s)
                """
                for tag in new_tags_list:
                    #Check if new tag is already in database
                    if tag.lower() not in cached_tags:
                        #New tag. Add it to database and keep its id
                        cursor.execute(add_tag_sql, (tag,))
                        cnx.commit()
                        tag_id = cursor.lastrowid
                    else:
                        #Tag exists. Get its id
                        tag_id = cached_tags[tag.lower()]

                    cursor.execute(add_datasettag_sql, (tag_id, dataset_id,) )
                    cnx.commit()

                # get the updated value...
                result['dataset'].extend(get_dataset(cursor, field2change, user.id, dataset_id))
                result['success'] = True

                cursor.close()
                cnx.close()

                print(json.dumps(result))


        elif field2change == "schematic":

            #get old schematic file
            old_file = get_old_schematic(cursor, user.id, dataset_id)

            #delete old schematic file if present
            if old_file is not None:
                result['old_schematic_file'] = 'Found and removed'
                os.remove( "{0}/{1}".format('..', old_file) )
            else:
                result['old_schematic_file'] = 'None found'

        	#set img's src path
            schematic_src_path = "{0}/{1}".format('../uploads/files', new_value)
            #print(json.dumps(new_value))

            #file is renamed with dataset uid. need to get the file extension
            match = re.search('\.(\w+)$', new_value)
            schematic_format = match.group(0)

            #set img's dest path
            schematic_dest_path = "{0}/{1}{2}".format('img/dataset_schematics', dataset_id, schematic_format)
            new_value = schematic_dest_path

            #update database entry
            change_dataset_info(cursor, field2change, user.id, dataset_id, new_value)
            cnx.commit()

            #move file from src path to dest path
            shutil.move(schematic_src_path, "../{0}".format(schematic_dest_path))

            # get the updated value...
            result['dataset'].extend(get_dataset(cursor, field2change, user.id, dataset_id))
            result['success'] = True

            cursor.close()
            cnx.close()

            print(json.dumps(result))

        else:
            result = { 'error':[] }

            error = "Not able to change dataset's information. Unexpected field was changed."
            result['error'] = error
            result['success'] = False

            print(json.dumps(result))

    else:
        result = { 'error':[] }

        error = "Not able to change dataset's information. You do not own this dataset."
        result['error'] = error
        result['owns_dataset'] = owns_dataset
        result['success'] = False

        print(json.dumps(result))


def check_dataset_ownership(cursor, current_user_id, dataset_id):
    qry = """
       SELECT d.id, d.owner_id
       FROM dataset d
       WHERE d.id = %s
    """
    cursor.execute(qry, (dataset_id,))

    # default: Assume user does not own dataset
    user_owns_dataset = False

    for row in cursor:
        # Change access if user owns the dataset
        if row[1] == current_user_id:
            user_owns_dataset = True

        # Return a statement that the user does not own the dataset (have permission)
        else:
        	user_owns_dataset = False

    return user_owns_dataset


def get_old_schematic(cursor, current_user_id, dataset_id):
    qry = """
       SELECT d.schematic_image
         FROM dataset d
              JOIN organism o ON d.organism_id=o.id
        WHERE d.id = %s
    """
    cursor.execute(qry, (dataset_id,))

    for row in cursor:
        old_filename = row[0]

    return old_filename


def change_dataset_info(cursor, field2change, current_user_id, dataset_id, new_value):

    if field2change == 'title':
        qry = """
            UPDATE dataset d
            JOIN guser g ON d.owner_id = g.id
            SET d.title = %s
            WHERE d.id = %s AND g.id = %s
        """
    if field2change == 'ldesc':
        qry = """
            UPDATE dataset d
            JOIN guser g ON d.owner_id = g.id
            SET d.ldesc = %s
            WHERE d.id = %s AND g.id = %s
        """
    if field2change == 'pubmed':
        qry = """
            UPDATE dataset d
            JOIN guser g ON d.owner_id = g.id
            SET	d.pubmed_id = %s
            WHERE d.id = %s AND g.id = %s
        """
    if field2change == 'geo':
        qry = """
            UPDATE dataset d
            JOIN guser g ON d.owner_id = g.id
            SET d.geo_id = %s
            WHERE d.id = %s AND g.id = %s
        """
    if field2change == 'schematic':
        qry = """
            UPDATE dataset d
            JOIN guser g ON d.owner_id = g.id
            SET	d.schematic_image = %s
            WHERE d.id = %s AND g.id = %s
        """

    cursor.execute(qry, (new_value, dataset_id, current_user_id,))


def get_dataset(cursor, field2change, current_user_id, dataset_id):

    if field2change == 'title':
        qry = """
           SELECT d.title
             FROM dataset d
                  JOIN organism o ON d.organism_id=o.id
            WHERE d.id = %s
        """
    if field2change == 'ldesc':
        qry = """
           SELECT d.ldesc
             FROM dataset d
                  JOIN organism o ON d.organism_id=o.id
            WHERE d.id = %s
        """
    if field2change == 'pubmed':
        qry = """
           SELECT d.pubmed_id
             FROM dataset d
                  JOIN organism o ON d.organism_id=o.id
            WHERE d.id = %s
        """
    if field2change == 'geo':
        qry = """
          SELECT d.geo_id
            FROM dataset d
                 JOIN organism o ON d.organism_id=o.id
          WHERE d.id = %s
        """
    if field2change == 'schematic':
        qry = """
           SELECT d.schematic_image
             FROM dataset d
                  JOIN organism o ON d.organism_id=o.id
            WHERE d.id = %s
        """
    if field2change == 'tag':
        qry = """
            SELECT IFNULL(GROUP_CONCAT(t.label), 'NULL')
              FROM tag t
                    LEFT JOIN dataset_tag d ON d.tag_id=IFNULL(t.id, 'NULL')
            WHERE d.dataset_id = %s;
        """

    cursor.execute(qry, (dataset_id,))
    dataset = list()

    for row in cursor:
        if field2change == 'pubmed' or field2change == 'geo' or field2change == 'tag':
            # displays 'Pubmed ID: null' or 'GEO ID: null' without this
            if row[0] is None:
                value = ''
            else:
                value = row[0]
            dataset.append({ 'updated_value': value })

        else:
            dataset.append({ 'updated_value': row[0] })

    return dataset


if __name__ == '__main__':
    main()
