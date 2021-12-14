#!/opt/bin/python3

"""
This can be used to change ownership of a single dataset or a 
profile and all the datasets within it.

If the --include_displays argument is passed, those displays owned
by the same original owner as the dataset are also updated.

Use with care.
"""
import argparse
import os
import sys

lib_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
sys.path.append(lib_path)
import geardb

def main():
    parser = argparse.ArgumentParser( description='Changes ownership of datasets or profiles')
    parser.add_argument('-no', '--new_owner_id', type=int, required=True, help='New owner guser.ID' )
    parser.add_argument('-id', '--include_displays', action='store_true', help='Changes configured displays also' )
    parser.add_argument('-d', '--dataset_id', type=str, required=False, help='Pass a single dataset.id' )
    parser.add_argument('-p', '--profile_id', type=str, required=False, help='Pass a profile/layout id' )
    args = parser.parse_args()

    if not args.dataset_id and not args.profile_id:
        print("ERROR: You must specify either --dataset_id or --profile_id (or both)")
    
    cnx = geardb.Connection()
    cursor = cnx.get_cursor()

    if args.dataset_id:
        reown_dataset(cursor, args.dataset_id, args.new_owner_id, args.include_displays)

    if args.profile_id:
        reown_profile(cursor, args.profile_id, args.new_owner_id, args.include_displays)

    cnx.commit()
    cursor.close()
    cnx.close()

def reown_dataset(cursor, dataset_id, user_id, include_displays):
    dataset = geardb.get_dataset_by_id(id=dataset_id)

    if dataset is None:
        raise Exception("ERROR: No dataset found of ID: {0}".format(dataset_id))
    
    original_owner = dataset.owner_id
    print("\nDataset:{0} Changing owner from {1} to {2}".format(dataset_id, original_owner, user_id))
    dataset.save_change(attribute='owner_id', value=user_id)

    ## check, are we doing curated displays too?
    if include_displays:
        qry = "SELECT id, user_id FROM dataset_display WHERE dataset_id = %s"
        cursor.execute(qry, (dataset_id,))

        display_ids_to_update = list()

        for (id, display_owner_id) in cursor:
            if display_owner_id == original_owner:
                display_ids_to_update.append(id)

        display_qry = "UPDATE dataset_display SET user_id = %s WHERE id = %s"
        pref_update_qry = """
             UPDATE dataset_preference 
             SET user_id = %s 
             WHERE user_id = %s
               AND dataset_id = %s
               AND display_id = %s
        """
        
        for display_id in display_ids_to_update:
            print("Dataset display:{0} Changing owner to {1}".format(display_id, user_id))
            cursor.execute(display_qry, (user_id, display_id))

            print("Dataset preference: Setting user_id to {0} where user_id was {1}, dataset_id {2} and display_id:{3}".format(user_id, original_owner, dataset_id, display_id))
            cursor.execute(pref_update_qry, (user_id, original_owner, dataset_id, display_id))
    

def reown_profile(cursor, profile_id, user_id, include_displays):
    layout = geardb.get_layout_by_id(profile_id)

    if layout is None:
        raise Exception("ERROR: No profile found of ID: {0}".format(profile_id))

    qry = "UPDATE layout SET user_id = %s WHERE id = %s"
    print("Layout:{0} Changing owner to {1}".format(profile_id, user_id))
    cursor.execute(qry, (user_id, profile_id))
    
    # now layout_members
    layout.get_members()

    print("Got {0} members in layout {1}".format(len(layout.members), profile_id))

    for lm in layout.members:
        reown_dataset(cursor, lm.dataset_id, user_id, include_displays)

if __name__ == '__main__':
    main()
    
