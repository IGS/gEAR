#!/opt/bin/python3

'''

To seek out duplicate shares:
select user_id, dataset_id, count(dataset_id), max(id)
  from dataset_shares
group by user_id, dataset_id
having count(dataset_id) > 1;

'''

import argparse
import mysql.connector
import configparser
import os
import re

def main():
    parser = argparse.ArgumentParser( description='')

    ## output file to be written
    parser.add_argument('-i', '--input_file', type=str, required=False, help='Path to an input file to be read' )
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read('gear.ini')

    try:
        cnx = mysql.connector.connect(user=config['database']['user'], password=config['database']['password'],
                                      host=config['database']['host'], database=config['database']['name'])

    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)

    cursor = cnx.cursor()

    # if false, no db changes are made and statements are printed instead
    debugging = True

    group_label = 'HRP'
    #group_label = 'NIDCD-Griffith'
    #group_label = 'Stone & Segil'
    layout_label = 'HRP'
    #layout_label = 'NIDCD-Griffith'
    #layout_label = 'Stone & Segil'

    # empty this list to do all members
    member_include_ids = []

    # get a list of people in this group
    member_ids = get_group_member_ids(cursor, group_label)

    # does this layout need to be added for each?  Comment this out if not
    #add_layouts(cursor, member_ids, layout_label)
    #cnx.commit()

    # for each of these, get their corresponding layout ID we want to add to
    member_layout_ids = get_layout_ids_by_member(cursor, member_ids, layout_label)

    # add each dataset to each member's layout
    #dataset_ids = ['bee735e5-d180-332c-7892-dd751dd76bb8'] # 'NIDCD-Griffith'
    # HRP
    dataset_ids = [
#         'e7c4268b-3c34-2eef-91ca-25f5c8c6b352',
#         '84878927-3f14-7831-a7f6-ec1e0402ed61',
#         '3ea4e4f7-c227-bc15-0af0-16d053251b3f',
#        'a7027971-1461-c8b6-299a-99f247139a50',
#        'a614bf59-1748-333b-d9ff-de635e58022e',
#         'ebd45945-342b-c865-25b3-82cc9eca2a19',
#         '5b8185ba-0a94-73cd-c6f0-fa706020def6',
#         '4c765b04-8a0e-e351-c08d-2ec312be2826',
#         'ede12ea0-9669-4af9-9786-3ba378b5b843',
#         '646ddd8d-183d-ad8a-fa3d-b206224fa7d8',
#         'liasdf99-1912-3oiu-h981-2bj0912g89h2',
#         'liasdf98-1129-s82n-slk1-ane128h21kad',
#         '2f4dc784-f581-6a43-0c51-0613b16c4930',
#         '19bb9a55-1815-0b37-2548-6e1e83f75cfd',
          'ab859cd1-2c4c-48a1-8ba0-0c0480e08f20'
    ]

    # Stone & Segil
    #dataset_ids = ['548a9372-084c-12f6-f18a-74c06ca6e547', 'ede12ea0-9669-4af9-9786-3ba378b5b843']

    for member_id in member_layout_ids:
        # are we filtering by ID?
        if len(member_include_ids) > 0:
            if member_id not in member_include_ids:
                continue

        print("Processing member ID: {0}".format(member_id))

        grid_position = get_layout_member_count(cursor, member_layout_ids[member_id]) + 1

        # Don't feel like figuring out math for proper start row and column, so let's start with a large row
        start_row = grid_position
        start_col = 1
        for dataset_id in dataset_ids:
            # TODO: skip this step if the user owns the dataset or if it has already been shared with them.
            qry = "INSERT INTO dataset_shares (dataset_id, user_id, is_allowed) VALUES (%s, %s, 1)"

            if debugging:
                print("{0} - ({1},{2})".format(qry, dataset_id, member_id))
            else:
                cursor.execute(qry, (dataset_id, member_id))

            # get the default single gene and multi gene display ID preferences for the dataset ID

            single_gene_preference = """
                SELECT display_id FROM dataset_preference
                WHERE dataset_id = %s AND is_multigene = 0
            """
            multi_gene_preference = """
                SELECT display_id FROM dataset_preference
                WHERE dataset_id = %s AND is_multigene = 1
            """

            cursor.execute(single_gene_preference, (dataset_id,))
            single_fetch = cursor.fetchone()

            cursor.execute(multi_gene_preference, (dataset_id,))
            multi_fetch = cursor.fetchone()

            # insert the new layout member record as separate single- and multi-gene display records
            new_layout_member_qry = """
                INSERT INTO layout_displays (layout_id, display_id, grid_position, start_col, grid_width, start_row, grid_height, math_preference)
                VALUES (%s, %s, %s, %s, %s, %s, %s, 'raw')
            """
            if single_fetch:
                single_gene_display_id = single_fetch[0]
                if debugging:
                    print("{0} - ({1},{2},{3},{4},{5})".format(new_layout_member_qry, member_layout_ids[member_id], single_gene_display_id, grid_position, 4, 0, 1))
                cursor.execute(new_layout_member_qry, (member_layout_ids[member_id], single_gene_display_id, grid_position, start_col, 4, start_row, 1))

            if multi_fetch:
                multi_gene_display_id = multi_fetch[0]
                if debugging:
                    print("{0} - ({1},{2},{3},{4},{5})".format(new_layout_member_qry, member_layout_ids[member_id], multi_gene_display_id, grid_position, 4, 0, 1))
                cursor.execute(new_layout_member_qry, (member_layout_ids[member_id], multi_gene_display_id, grid_position, start_col, 4, start_row, 1))

            start_col += 4
            if start_col > 12:
                start_col = 1
                start_row += 1

            grid_position += 1

    cnx.commit()
    cursor.close()
    cnx.close()

def add_layouts(curs, ids, label):
    for user_id in ids:
        add_sql = """
           INSERT INTO layout (user_id, label, is_current)
           VALUES (%s, %s, 0)
        """
        if debugging:
            print("Here I'd add layout {0}".format(label))
        else:
            curs.execute(add_sql, (user_id, label))


def get_group_member_ids(curs, label):
    qry = """SELECT gm.user_id
               FROM ggroup_members gm
                    JOIN ggroup g ON g.id=gm.group_id
              WHERE g.label = %s
    """
    member_ids = list()

    curs.execute(qry, (label,))

    for (user_id,) in curs:
        member_ids.append(user_id)

    return member_ids

def get_layout_ids_by_member(curs, member_ids, layout_label):
    member_layouts = dict()

    for member_id in member_ids:
        qry = """
            SELECT id FROM layout WHERE user_id = %s AND label = %s
        """
        curs.execute(qry, (member_id, layout_label))

        layout_found = False
        for (layout_id,) in curs:
            layout_found = True
            member_layouts[member_id] = layout_id

        if layout_found is False:
            raise Exception("Error, No layout called {0} found for user {1}".format(layout_label, member_id))

    return member_layouts

def get_layout_member_count(curs, layout_id):
    qry = "SELECT count(*) FROM layout_displays WHERE layout_id = %s"

    # grid position will be higher than it should be because layout members are not split
    #  into single and multi gene displays.  This will not affect the order.

    curs.execute(qry, (layout_id,))

    for (member_count,) in curs:
        return member_count

if __name__ == '__main__':
    main()







