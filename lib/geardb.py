import datetime
import json
import os
import re
import sys
import uuid

from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Optional
from json import JSONEncoder

gear_lib_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(gear_lib_path)
import gear.db

# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]

from gear.serverconfig import ServerConfig
this.servercfg = ServerConfig().parse()

# This is where things specific to dynamic analyses will be stored, such as intermediate
#  H5AD files, images, etc.
this.analysis_base_dir = '/tmp'

# Overrides the json module so JSONEncoder.default() automatically checks for to_json()
#  in any class to be directly serializable.
#  Ref: https://stackoverflow.com/a/38764817/1368079
def _default(self, obj):
    try:
        return getattr(obj.__class__, "_serialize_json", _default.default)(obj)
    except:
        return str(obj)

_default.default = JSONEncoder().default
JSONEncoder.default = _default

def _read_site_domain_config():
    """Convert site domain preferences into a dictionary."""
    this_dir = os.path.dirname(os.path.abspath(__file__))
    with open ("{0}/../www/site_domain_prefs.json".format(this_dir)) as json_file:
        return json.loads(json_file.read())

def _read_domain_url():
    json_conf = _read_site_domain_config()
    # Can build a longer URL off of this one
    return json_conf['domain_url']

def _read_domain_label():
    json_conf = _read_site_domain_config()
    return json_conf['domain_label']

def _read_domain_links_out():
    json_conf = _read_site_domain_config()
    return json_conf['links_out']

def _read_domain_short_label():
    json_conf = _read_site_domain_config()
    return json_conf['domain_short_display_label']



# For those functional differences we have depending on the site domain
#  Current values are gear, nemo.  Value kept in www/site_domain_prefs.json
# TODO: each of these, as currently implemented, causes file I/O and shouldn't.
this.domain_url = _read_domain_url()
this.domain_label = _read_domain_label()
this.domain_short_label = _read_domain_short_label()
this.links_out = _read_domain_links_out()

def check_verification_code(long_form, short_form):
    """
    To prevent fake users from creating accounts, we now employ e-mail verification.
    We don't want to store a temp verification code in the database, so we generate
    a full UUID in javascript and then hash it to a short form for user's verification
    e-mail.  This function checks if the long form matches the short form.
    """
    if get_verification_code_short_form(long_form) == short_form:
        return True
    else:
        return False

def get_verification_code_short_form(long_form):
    """
    Get the hashed/short form of the full verification code. If you want to ensure
    anyone who can see the repo in GitHub can't script account creation still, you
    can just modify the string this function returns for your own site.
    """
    return ''.join([x[0] for x in long_form.split('-')])

def get_analysis(analysis, dataset_id, session_id, is_spatial=False):
    """Return analysis object based on various factors."""
    # If an analysis is posted we want to read from its h5ad
    if analysis:
        user = get_user_from_session_id(session_id)
        user_id = None
        if user:
            user_id = user.id

        if is_spatial:
            ana = SpatialAnalysis(id=analysis['id'], dataset_id=dataset_id,
                                session_id=session_id, user_id=user_id)
        else:
            ana = Analysis(id=analysis['id'], dataset_id=dataset_id,
                                session_id=session_id, user_id=user_id)

        if 'type' in analysis:
            ana.type = analysis['type']
        else:
            ana.discover_type()

        # Check that the h5ad file exists
        if not os.path.exists(ana.dataset_path()):
            raise FileNotFoundError("No h5 file found for the passed in analysis {}".format(ana.dataset_path()))

    else:

        ds = Dataset(id=dataset_id, has_h5ad=1)
        filetype = "h5"

        # Ensure the zarr file is retrieved instead of the h5
        if is_spatial:
            ds.dtype = 'spatial'
            filetype = "zarr"
            ds.has_h5ad = 0 # does not affect get_file_path but sanity-checking

        h5_path = ds.get_file_path()

        # Let's not fail if the file isn't there
        if not os.path.exists(h5_path):
            raise FileNotFoundError("No {} file found for this dataset {}".format(filetype, h5_path))
        ana = Analysis(type='primary', dataset_id=dataset_id)
    return ana

def get_dataset_by_id(id=None, include_shape=None):
    """
    Given a dataset ID string this returns a Dataset object with all attributes
    populated which come directly from the table.  Secondary things, such as tags,
    must be loaded separately.

    Returns None if no dataset of that ID is found.
    """

    conn = Connection()
    cursor = conn.get_cursor()

    qry = """
         SELECT id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, is_downloadable, ldesc, date_added,
                dtype, schematic_image, share_id, math_default, marked_for_removal, load_status,
                has_h5ad
           FROM dataset
          WHERE id = %s
    """
    cursor.execute(qry, (id, ))
    dataset = None

    for (id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, is_downloadable, ldesc, date_added,
         dtype, schematic_image, share_id, math_default, marked_for_removal, load_status,
         has_h5ad) in cursor:
        dataset = Dataset(id=id, owner_id=owner_id, title=title, organism_id=organism_id,
                          pubmed_id=pubmed_id, geo_id=geo_id, is_public=is_public, is_downloadable=is_downloadable, ldesc=ldesc,
                          date_added=date_added, dtype=dtype, schematic_image=schematic_image,
                          share_id=share_id, math_default=math_default,
                          marked_for_removal=marked_for_removal, load_status=load_status,
                          has_h5ad=has_h5ad)

    if include_shape == '1':
        dataset.get_shape()

    cursor.close()
    conn.close()
    return dataset

def get_dataset_by_share_id(share_id=None, include_shape=None):
    """
    Given a dataset share ID string this returns a Dataset object with all attributes
    populated which come directly from the table.  Secondary things, such as tags,
    must be loaded separately.

    Returns None if no dataset of that shareID is found.
    """

    conn = Connection()
    cursor = conn.get_cursor()

    qry = """
         SELECT id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, is_downloadable, ldesc, date_added,
                dtype, schematic_image, share_id, math_default, marked_for_removal, load_status,
                has_h5ad
           FROM dataset
          WHERE share_id = %s
    """
    cursor.execute(qry, (share_id, ))
    dataset = None

    for (id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, is_downloadable, ldesc, date_added,
         dtype, schematic_image, share_id, math_default, marked_for_removal, load_status,
         has_h5ad) in cursor:
        dataset = Dataset(id=id, owner_id=owner_id, title=title, organism_id=organism_id,
                          pubmed_id=pubmed_id, geo_id=geo_id, is_public=is_public, is_downloadable=is_downloadable, ldesc=ldesc,
                          date_added=date_added, dtype=dtype, schematic_image=schematic_image,
                          share_id=share_id, math_default=math_default,
                          marked_for_removal=marked_for_removal, load_status=load_status,
                          has_h5ad=has_h5ad)

    if include_shape == '1':
        dataset.get_shape()

    cursor.close()
    conn.close()
    return dataset

def get_dataset_by_title(title=None, include_shape=None):
    """
    Given a dataset title string this returns a Dataset object with all attributes
    populated which come directly from the table.  Secondary things, such as tags,
    must be loaded separately.

    Returns None if no dataset of that title is found.

    Throws an exception if more than one dataset with that same title is found
    """

    conn = Connection()
    cursor = conn.get_cursor()

    qry = """
         SELECT id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, is_downloadable, ldesc, date_added,
                dtype, schematic_image, share_id, math_default, marked_for_removal, load_status,
                has_h5ad
           FROM dataset
          WHERE title = %s
            AND marked_for_removal = 0
    """
    cursor.execute(qry, (title, ))
    dataset = None

    found = 0

    for (id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, is_downloadable, ldesc, date_added,
         dtype, schematic_image, share_id, math_default, marked_for_removal, load_status,
         has_h5ad) in cursor:
        dataset = Dataset(id=id, owner_id=owner_id, title=title, organism_id=organism_id,
                          pubmed_id=pubmed_id, geo_id=geo_id, is_public=is_public, is_downloadable=is_downloadable, ldesc=ldesc,
                          date_added=date_added, dtype=dtype, schematic_image=schematic_image,
                          share_id=share_id, math_default=math_default,
                          marked_for_removal=marked_for_removal, load_status=load_status,
                          has_h5ad=has_h5ad)
        found += 1

    cursor.close()
    conn.close()

    if found == 1:
        if include_shape == '1':
            dataset.get_shape()

        return dataset
    elif found == 0:
        return None
    raise Exception("Error: More than one dataset found with the same title: {0}".format(dataset.title))

def get_dataset_collection(ids=None):
    return 1

def get_dataset_count():
    conn = Connection()
    cursor = conn.get_cursor()

    qry = "SELECT count(id) FROM dataset WHERE marked_for_removal = 0"
    cursor.execute(qry)

    for (c) in cursor:
        cursor.close()
        conn.close()
        return c[0]

def get_dataset_id_from_share_id(share_id):
    """
    Given a share_id passed this returns a numeric dataset.id if a corresponding one
    is found. Otherwise, None is returned
    """
    conn = Connection()
    cursor = conn.get_cursor()
    dataset_id = None

    qry = "SELECT id FROM dataset WHERE share_id = %s"
    cursor.execute(qry, (share_id,))

    for row in cursor:
        dataset_id = row[0]

    cursor.close()
    conn.close()

    return dataset_id


def get_gene_by_id(gene_id):
    """
    Given a gene_id passed this returns a Gene object with all attributes populated. Returns
    None if no gene is found with that ID.
    """
    conn = Connection()
    cursor = conn.get_cursor()

    qry = """
          SELECT g.id, g.ensembl_id, g.genbank_acc, g.organism_id, gs.label, g.product, g.biotype
            FROM gene g
                 JOIN gene_symbol gs ON gs.gene_id=g.id
           WHERE gs.is_primary = 1
             AND g.id = %s
    """
    cursor.execute(qry, (gene_id,))
    gene = None

    for (id, ensembl_id, genbank_acc, organism_id, gene_symbol, product, biotype) in cursor:
        gene = Gene(id=id, ensembl_id=ensembl_id, genbank_acc=genbank_acc, organism_id=organism_id,
                    gene_symbol=gene_symbol, product=product, biotype=biotype)
        break

    cursor.close()
    conn.close()
    return gene

def get_gene_cart_by_id(gc_id):
    """
    Given a gene_cart_id passed this returns a GeneCart object with all attributes populated. Returns
    None if no cart is found with that ID.
    """
    conn = Connection()
    cursor = conn.get_cursor()

    qry = """
          SELECT id, user_id, organism_id, gctype, label, ldesc, share_id, is_public, date_added
            FROM gene_cart
           WHERE id = %s
    """
    cursor.execute(qry, (gc_id,))
    gc = None

    for (id, user_id, organism_id, gctype, label, ldesc, share_id, is_public, date_added) in cursor:
        gc = GeneCart(id=id, user_id=user_id, organism_id=organism_id, gctype=gctype, label=label, ldesc=ldesc,
                      share_id=share_id, is_public=is_public, date_added=date_added)
        break

    cursor.close()
    conn.close()
    return gc

def get_gene_cart_by_share_id(share_id):
    """
    Given a gene_cart_id passed this returns a GeneCart object with all attributes populated. Returns
    None if no cart is found with that ID.
    """
    conn = Connection()
    cursor = conn.get_cursor()

    qry = """
          SELECT id, user_id, organism_id, gctype, label, ldesc, share_id, is_public, date_added
            FROM gene_cart
           WHERE share_id = %s
    """
    cursor.execute(qry, (share_id,))
    gc = None

    for (id, user_id, organism_id, gctype, label, ldesc, share_id, is_public, date_added) in cursor:
        gc = GeneCart(id=id, user_id=user_id, organism_id=organism_id, gctype=gctype, label=label, ldesc=ldesc,
                      share_id=share_id, is_public=is_public, date_added=date_added)
        break

    cursor.close()
    conn.close()
    return gc

def get_layout_by_id(layout_id):
    """
    Given a passed layout_id returns a Layout object with all attributes
    populated.  Returns None if no layout is found with that ID.
    """
    conn = Connection()
    cursor = conn.get_cursor()
    layout = None

    qry = """
          SELECT id, user_id, label, is_current, is_domain, share_id, is_public
          FROM layout
          WHERE id = %s
    """
    cursor.execute(qry, (layout_id,))

    for (id, user_id, label, is_current, is_domain, share_id, is_public) in cursor:
        layout = Layout(id=id, user_id=user_id, is_domain=is_domain,
                        label=label, is_current=is_current,
                        share_id=share_id, is_public=is_public)
        break

    cursor.close()
    conn.close()
    return layout

def get_layout_by_share_id(layout_share_id):
    """
    Given a passed layout_share_id returns a Layout object with all attributes
    populated.  Returns None if no layout is found with that ID.
    """
    conn = Connection()
    cursor = conn.get_cursor()
    layout = None

    qry = """
          SELECT id, user_id, label, is_current, is_domain, share_id, is_public
          FROM layout
          WHERE share_id = %s
    """
    cursor.execute(qry, (layout_share_id,))

    for (id, user_id, label, is_current, is_domain, share_id, is_public) in cursor:
        layout = Layout(id=id, user_id=user_id, is_domain=is_domain,
                        label=label, is_current=is_current,
                        share_id=share_id, is_public=is_public)
        break

    cursor.close()
    conn.close()
    return layout

def get_user_count():
    conn = Connection()
    cursor = conn.get_cursor()

    qry = "SELECT count(id) FROM guser"
    cursor.execute(qry)

    for (c) in cursor:
        cursor.close()
        conn.close()
        return c[0]

def get_user_by_id(user_id) -> "User | None":
    """
       Given a user_id string this returns a User object with
       all attributes populated.  Returns None if no user with
       that ID is found.
    """
    conn = Connection()
    cursor = conn.get_cursor()

    qry = """
          SELECT g.id, g.user_name, g.email, g.institution, g.pass, g.updates_wanted,
                 g.is_admin, g.default_org_id, g.is_curator, g.help_id, g.colorblind_mode,
                 l.share_id
            FROM guser g
                 LEFT JOIN layout l ON g.layout_id=l.id
           WHERE g.id = %s
    """
    cursor.execute(qry, (user_id, ) )

    user = None
    row = cursor.fetchone()
    if row:
        (id, user_name, email, institution, password, updates_wanted, is_admin,
         default_org_id, is_curator, help_id, colorblind_mode, layout_share_id) = row
        user = User(id=id, user_name=user_name, email=email, institution=institution,
                    password=password, updates_wanted=updates_wanted, is_admin=is_admin,
                    default_org_id=default_org_id, is_curator=is_curator, help_id=help_id,
                    colorblind_mode=colorblind_mode, layout_share_id=layout_share_id)

    cursor.close()
    conn.close()
    return user

def get_user_from_session_id(session_id: str | None) -> "User | None":
    """
    Retrieve a User object based on the provided session ID.

    Args:
        session_id (str | None): The session ID associated with the user session.

    Returns:
        User | None: A User object if a matching session is found; otherwise, None.

    Queries the database to find a user associated with the given session ID by joining
    the `guser`, `user_session`, and `layout` tables. Returns a User object populated
    with user details if a valid session is found, or None if not.

    """

    conn = Connection()
    cursor = conn.get_cursor()

    if not session_id:
        cursor.close()
        conn.close()
        return None

    qry = """
          SELECT g.id, g.user_name, g.email, g.institution, g.pass, g.updates_wanted,
                 g.is_admin, g.default_org_id, g.is_curator, g.help_id, g.colorblind_mode,
                 l.share_id
            FROM guser g
                 JOIN user_session us ON g.id=us.user_id
                 LEFT JOIN layout l ON g.layout_id=l.id
           WHERE us.session_id = %s
    """
    cursor.execute(qry, (session_id, ) )

    user = None
    for (id, user_name, email, institution, password, updates_wanted, is_admin, default_org_id,
         is_curator, help_id, colorblind_mode, layout_share_id) in cursor:
        user = User(id=id, user_name=user_name, email=email, institution=institution,
                    password=password, updates_wanted=updates_wanted, is_admin=is_admin,
                    default_org_id=default_org_id, is_curator=is_curator, help_id=help_id,
                    colorblind_mode=colorblind_mode, layout_share_id=layout_share_id)
        break

    cursor.close()
    conn.close()
    return user

def get_display_by_id(display_id):
    """Given user id, return all datasets representations from dataset display table."""
    conn = Connection()
    cursor = conn.get_cursor()
    qry = "SELECT * from dataset_display where id = %s"
    cursor.execute(qry, (display_id,))

    try:
        (id, dataset_id, user_id, label, plot_type, plotly_config) = cursor.fetchone()
        return dict(
            id=id, dataset_id=dataset_id, user_id=user_id, label=label,
            plot_type=plot_type, plotly_config=json.loads(plotly_config)
        )
    except:
        return None
    finally:
        cursor.close()
        conn.close()

def get_default_display(user_id, dataset_id, is_multigene=0):
    """Return user's display preference for given dataset."""
    conn = Connection()
    cursor = conn.get_cursor()
    qry = """
        SELECT display_id FROM dataset_preference
        where user_id = %s and dataset_id = %s and is_multigene = %s
    """
    cursor.execute(qry, (user_id, dataset_id, is_multigene))
    try:
        (default_display_id,) = cursor.fetchone()
        return default_display_id
    except:
        # User has no display preference for this dataset
        return None
    finally:
        cursor.close()
        conn.close()

def get_displays_by_user_id(user_id, dataset_id):
    """Given user id, return all datasets representations from dataset display table."""
    conn = Connection()
    cursor = conn.get_cursor()
    qry = "SELECT * from dataset_display where user_id = %s and dataset_id = %s"
    cursor.execute(qry, (user_id, dataset_id))
    displays = [
        dict(
            id=id, dataset_id=dataset_id, user_id=user_id, label=label,
            plot_type=plot_type, plotly_config=json.loads(plotly_config)
        )
        for (id, dataset_id, user_id, label, plot_type, plotly_config) in cursor
     ]
    cursor.close()
    conn.close()
    return displays

def get_user_id_from_session_id(session_id):
    """Return the user id for the given session."""
    conn = Connection()
    cursor = conn.get_cursor()
    qry = ( "SELECT user_id FROM user_session WHERE session_id = %s" )
    cursor.execute(qry, (session_id, ) )
    user_id = None

    for (uid,) in cursor:
        user_id = uid
    cursor.close()
    conn.close()
    return user_id


def get_gene_by_gene_symbol(gene_symbol, dataset_id):
    qry_org_id = "SELECT organism_id from dataset where id = %s"

    conn = Connection()
    cursor = conn.get_cursor()
    cursor.execute(qry_org_id, (dataset_id,))

    (org_id,) = cursor.fetchone()

    qry_gene_location = """
            SELECT id, ensembl_id, ensembl_version, ensembl_release, genbank_acc, organism_id, molecule, start, stop, gene_symbol, product, biotype
            FROM gene
            WHERE gene_symbol = %s and organism_id = %s order by ensembl_release DESC
        """

    cursor.execute(qry_gene_location, (gene_symbol, org_id,))

    (id, ensembl_id, ensembl_version, ensembl_release, genbank_acc, organism_id, molecule, start, stop, gene_symbol, product, biotype,) = cursor.fetchone()

    cursor.close()
    conn.close()

    gene = Gene(id=id, ensembl_id=ensembl_id, ensembl_version=ensembl_version,
                ensembl_release=ensembl_release, genbank_acc=genbank_acc,
                organism_id=organism_id, molecule=molecule, start=start,
                stop=stop, gene_symbol=gene_symbol, product=product,
                biotype=biotype)

    return gene


class Analysis:
    """
    When printed directly the JSON representation of the analysis is returned.

    Path conventions:

    PRIMARY
    - These are the direct h5 files created by the user when they upload a dataset.
    -------
    www/datasets/$dataset_id.h5ad

    PUBLIC
    - These can only be made by owner or gear curators and are publicly available for everyone
      to see/copy.
    ------
    www/analyses/by_dataset/$dataset_id/$analysis_id/$dataset_id.h5ad

    USER_SAVED
    - Created by users from their datasets or any other public ones. Visible only to the user.
    ----------
    www/analyses/by_user/$user_id/$dataset_id/$analysis_id/$dataset_id.h5ad

    USER_UNSAVED
    - These are created automatically by the interface any time a user does an analysis step,
      saving progress.
    ------------
    /tmp/$session/$dataset_id/$analysis_id/$dataset_id.h5ad
    /tmp/e385305c-4387-433e-8b62-4bcf7c30ac52/ab859cd1-2c4c-48a1-8ba0-0c0480e08f20

    If a user selects a PRIMARY or PUBLIC analysis and makes modifications, it should first
    be copied to USER_UNSAVED, issued a new analysis_id, then changes made.

    Selecting a USER_SAVED or USER_UNSAVED should allow modifications directly.
    """

    def __init__(self, id=None, dataset_id=None, user_id=None, session_id=None, label=None, type=None, vetting=None):
        self.id = id
        self.dataset_id = dataset_id
        self.label = label
        self.session_id = session_id
        self.user_id = user_id

        # types are 'primary', 'public', 'user_saved', 'user_unsaved'
        self.type = type if type is not None else 'primary'

        # vettings are None, 'owner', 'gear', or 'community'
        self.vetting = vetting

        # if user ID wasn't set but the session was, do the lookup
        if user_id is None and session_id:
            self.user_id = get_user_id_from_session_id(session_id)

        if self.type not in ['primary', 'public', 'user_saved', 'user_unsaved']:
            raise Exception("ERROR: Invalid type '{0}' for Analysis instance".format(self.type))

        if self.vetting not in [None, 'owner', 'gear', 'community']:
            raise Exception("ERROR: Invalid vetting '{0}' for Analysis instance".format(self.vetting))

    def __repr__(self):
        pipeline_file = self.settings_path()
        json_data = json.loads(open(pipeline_file).read())
        # change "user_session_id" to "session_id" for consistency
        json_data['session_id'] = json_data.pop('user_session_id')
        return json.dumps(json_data, indent=4)

    def _serialize_json(self):
        # Called when json modules attempts to serialize
        return self.__dict__

    def base_path(self):
        """
        Returns the base directory path for an analysis, based on its actual type.  This allows for the support
        of parallel types 'primary', 'public', 'user_saved', 'user_unsaved'
        """
        this_dir = os.path.dirname(os.path.abspath(__file__))

        if self.type == 'primary':
            return "{0}/../www/datasets".format(this_dir)

        else:
            # all other types require analysis ID to be set
            if self.id is None:
                raise Exception("ERROR: base_path() called on Analysis object with no id attribute set.")

            if self.type == 'public':
                # ./$dataset_id/$analysis_id
                return "{0}/../www/analyses/by_dataset/{1}/{2}".format(this_dir, self.dataset_id, self.id)

            elif self.type == 'user_saved':
                if self.user_id is None:
                    raise Exception("ERROR: base_path() called on Analysis object with no user_id attribute set. Probably not logged in.")

                # ./$user_id/$dataset_id/$analysis_id/$dataset_id.h5ad
                return "{0}/../www/analyses/by_user/{1}/{2}/{3}".format(this_dir, self.user_id, self.dataset_id, self.id)

            elif self.type == 'user_unsaved':
                if self.session_id is None:
                    raise Exception("ERROR: base_path() called on Analysis object with no session_id attribute set. Probably not logged in.")

                # /tmp/$session/$dataset_id/$analysis_id/$dataset_id.h5ad
                return "/tmp/{0}/{1}/{2}".format(self.session_id, self.dataset_id, self.id)

    def dataset_path(self):
        return "{0}/{1}.h5ad".format(self.base_path(), self.dataset_id)

    def discover_vetting(self, current_user_id=None):
        """
        This describes the public attribution of the analysis.  Making it a derived value via this
        method rather an than explicitly stored one in the database or JSON layer so that changes
        made to the dataset ownership won't require updating of this as well.  This method will
        just return the correct value.

        Returns one of the values 'owner', 'gear' or 'community'

        If the owner is also a curator it gives priority to the curator status (gear)
        """

        # Analysis.user_id must be knownor we can't do this
        if self.user_id is None:
            return Exception("ERROR: Attempted to call Analysis.discover_vetting() without an owner assigned to the analysis")

        current_user = get_user_by_id(current_user_id)

        # if the user who created it is a gear curator, return that
        if current_user.is_curator:
            self.vetting = 'gear'
            return 'gear'
        # Are the current user and dataset owner the same?
        elif self.user_id == current_user_id:
            self.vetting = 'owner'
            return 'owner'
        else:
            self.vetting = 'community'
            return 'community'


    def discover_type(self):
        """
        Given an analysis ID it's technically possible to scan the directory hierarchies and
        find the type.

        Requires these attributes to be set:
        - dataset_id
        - user_id
        - session_id (if type is 'user_unsaved')

        Returns the discovered type AND sets it as self.type.

        Logic:
        1. If the current user is passed:
           1.1 -
        2. If the current user is not passed

        Returns None if not found
        """
        this_dir = os.path.dirname(os.path.abspath(__file__))

        # if the analysis ID and dataset ID are the same, it's a primary analysis
        if self.id == self.dataset_id:
            self.type = 'primary'
            return 'primary'

        # check user_saved
        test_path = "{0}/../www/analyses/by_user/{1}/{2}/{3}/{2}.h5ad".format(this_dir, self.user_id, self.dataset_id, self.id)
        if os.path.exists(test_path):
            self.type = 'user_saved'
            return 'user_saved'

        # check user_unsaved
        test_path = "/tmp/{0}/{1}/{2}/{1}.h5ad".format(self.session_id, self.dataset_id, self.id)
        if os.path.exists(test_path):
            self.type = 'user_unsaved'
            return 'user_unsaved'

        # Check for public first
        test_path = "{0}/../www/analyses/by_dataset/{1}/{2}/{1}.h5ad".format(this_dir, self.dataset_id, self.id)
        if os.path.exists(test_path):
            self.type = 'public'
            return 'public'

        # Didn't find it if we got this far
        return None

    @classmethod
    def from_json(cls, jsn):
        """
        Returns an Analyis object from a JSON object with the same attributes
        """

        # Dataset ID is now saved only in the dataset object, but some analyses may have it as a top-level attribute
        try:
            dataset_id = jsn["dataset"]["id"]
        except:
            dataset_id = jsn["dataset_id"]

        try:
            session_id = jsn["user_session_id"]
        except:
            session_id = jsn["analysis_session_id"]

        ana = Analysis(
            id=jsn['id'], dataset_id=dataset_id, label=jsn['label'], session_id=session_id,
            user_id=None, type=jsn['type']
        )

        ## get the rest of the properties
        for k in jsn:
            if not hasattr(ana, k):
                # some were manually named, skip them
                if k not in ['user_session_id', 'analysis_session_id']:
                    setattr(ana, k, jsn[k])

        return ana

    def get_adata(self, backed=False, force_sparse=False):
        """
        Returns the anndata object for the current analysis.
        """
        # This goes against PEP8, but putting the import here should speed up the
        #  common case of most of the rest of this module where scanpy isn't needed
        import scanpy as sc

        kwargs = {}

        if backed:
            kwargs['backed'] = 'r'
        if force_sparse:
            kwargs['as_sparse'] = "raw.X"

        return sc.read_h5ad(self.dataset_path(), **kwargs)


    def marker_gene_json_path(self):
        return "{0}/{1}.marker_gene_table.json".format(self.base_path(), self.dataset_id)

    def parent_path_by_type(self, type=None):
        """
        Returns the base directory base path for an analysis, based on a hypothetical type.  This allows for the support
        of parallel types 'primary', 'public', 'user_saved', 'user_unsaved'.

        This is generally useful if you want to find where an analysis directory would be if it existed, such as
        when searching for lists of analyses.
        """
        this_dir = os.path.dirname(os.path.abspath(__file__))

        if type == 'primary':
            return "{0}/../www/datasets".format(this_dir)

        else:
            if type == 'public':
                return "{0}/../www/analyses/by_dataset/{1}".format(this_dir, self.dataset_id)

            elif type == 'user_saved':
                if self.user_id is None:
                    raise Exception("ERROR: _parent_path_by_type() called on Analysis object with no user_id attribute set.")

                return "{0}/../www/analyses/by_user/{1}/{2}".format(this_dir, self.user_id, self.dataset_id)

            elif type == 'user_unsaved':
                if self.session_id is None:
                    raise Exception("ERROR: _parent_path_by_type() called on Analysis object with no session_id attribute set.")

                # /tmp/$session/$dataset_id/$analysis_id/$dataset_id.h5ad
                return "/tmp/{0}/{1}".format(self.session_id, self.dataset_id)

    def settings_path(self):
        return "{0}/{1}.pipeline.json".format(self.base_path(), self.dataset_id)

class SpatialAnalysis(Analysis):
    """
    Spatial-based analysis object.  Inherits from Analysis and adds spatial-specific methods.
    """

    def __init__(self, id=None, dataset_id=None, user_id=None, session_id=None, label=None, type=None, vetting=None):
        super().__init__(id=id, dataset_id=dataset_id, user_id=user_id, session_id=session_id, label=label, type=type, vetting=vetting)

    def dataset_path(self):
        return "{0}/spatial/{1}.zarr".format(self.base_path(), self.dataset_id)

    def determine_platform(self, sdata):
        try:
            platform = sdata.tables["table"].uns["platform"]
            return platform
        except KeyError:
            raise ValueError("No platform information found in the dataset")

    def discover_type(self):
        return super().discover_type()

    def get_sdata(self):
        import spatialdata as sd

        if not os.path.exists(self.dataset_path()):
            raise ValueError(f"Dataset {self.dataset_id} not found")
        return sd.read_zarr(self.dataset_path())

    def settings_path(self):
        base_path = f"{self.base_path()}/spatial"
        return "{0}/{1}.pipeline.json".format(base_path, self.dataset_id)

class AnalysisCollection:
    def __init__(self, public=None, user_saved=None, user_unsaved=None):
        self.public = [] if public is None else public
        self.user_saved = [] if user_saved is None else user_saved
        self.user_unsaved = [] if user_unsaved is None else user_unsaved

    def _scan_analysis_directory(self, ana, atype):
        """
        Searches a directory to find stored analysis files.  Assumes the dir passed is a parent directory which
        can contain more than one analysis directory.  Looking essentially for this:

        $dir/*/*.pipeline.json

        Returns a list of JSON objects, one for each analysis
        """
        analyses = list()
        dir = ana.parent_path_by_type(type=atype)

        if os.path.exists(dir):
            if atype == 'primary':
                ana.type = 'primary'
                json_path = ana.settings_path()

                if os.path.exists(json_path):
                    json_obj = json.loads(open(json_path).read())
                    analyses.append(Analysis.from_json(json_obj))
            else:
                for thing in os.listdir(dir):
                    dir_path = "{0}/{1}".format(dir, thing)

                    if os.path.isdir(dir_path):
                        for pipeline_file in (f for f in os.listdir(dir_path) if f.endswith('.pipeline.json')):
                            json_path = "{0}/{1}".format(dir_path, pipeline_file)
                            json_obj = json.loads(open(json_path, encoding="utf-8").read())
                            analyses.append(Analysis.from_json(json_obj))

        return analyses

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        # Called when json modules attempts to serialize
        return self.__dict__

    def get_all_by_dataset_id(self, user_id=None, session_id=None, dataset_id=None):
        """
        Gets all possible analyses for a given dataset, including primary, public, user-saved and
        user-unsaved analyses.
        """
        # clear any existing ones first
        self.__init__()

        ## Create hypothetical analysis to get paths
        ana = Analysis(dataset_id=dataset_id, user_id=user_id, session_id=session_id)

        ## Each of these is a list of JSON objects
        self.primary = self._scan_analysis_directory(ana, 'primary')
        self.public = self._scan_analysis_directory(ana, 'public')
        if user_id:
            self.user_saved = self._scan_analysis_directory(ana, 'user_saved')
        if session_id:
            self.user_unsaved = self._scan_analysis_directory(ana, 'user_unsaved')

class Connection:
    def __init__(self):
        self.mysql_cnx = gear.db.MySQLDB().connect()

    def commit(self):
        self.mysql_cnx.commit()

    def close(self):
        if self.mysql_cnx.is_connected():
            self.mysql_cnx.close()

    def get_cursor(self, use_dict=False):
        if use_dict == True:
            return self.mysql_cnx.cursor(dictionary=True)
        else:
            return self.mysql_cnx.cursor()

    def __del__(self):
        self.close()

class Organism:
    def __init__(self, id=None, label=None, genus=None, species=None, strain=None, taxon_id=None):
        self.id = id
        self.label = label
        self.genus = genus
        self.species = species
        self.strain = strain
        self.taxon_id = taxon_id

    def __repr__(self):
        return json.dumps(self.__dict__)

@dataclass
class OrganismCollection:
    organisms: List[Organism] = field(default_factory=list)

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        # Called when json modules attempts to serialize
        return self.__dict__

    def get_all(self):
        """
        Gets the full collection of organisms within the database, often used to generate
        select boxes on the interface.

        Populates the self.organisms list and returns the actual list of organisms
        """
        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              SELECT id, label, genus, species, strain, taxon_id
                FROM organism
            ORDER BY label
        """
        cursor.execute(qry)

        for row in cursor:
            org = Organism(id=row[0], label=row[1], genus=row[2], species=row[3],
                           strain=row[4], taxon_id=row[5]
            )
            self.organisms.append(org)

        cursor.close()
        conn.close()

        return self.organisms

class Layout:
    def __init__(self, id=None, user_id=None, is_domain=None, label=None,
                 is_current=None, share_id=None, members=[], folder_id=None,
                 folder_parent_id=None, folder_label=None, is_public=None):
        self.id = id
        self.user_id = user_id
        self.label = label
        self.is_current = is_current
        self.is_domain = is_domain
        self.is_public = is_public
        self.share_id = share_id

        # The are derived, populated by LayoutCollection methods
        self.folder_id = folder_id
        self.folder_parent_id = folder_parent_id
        self.folder_label = folder_label

        self.dataset_count = 0

        self.singlegene_members = list()
        self.multigene_members = list()

        # This should be a list of LayoutDisplays objects
        self.members = list()
        if members:
            self.members = members

        # handle defaults
        # TODO: If is_current = 1 we really need to reset all the other layouts by this user
        if not is_current:
            self.is_current = 0

        if not is_domain:
            self.is_domain = 0

        if not share_id:
            self.share_id = str(uuid.uuid4()).split('-')[0]

    def __repr__(self):
        return json.dumps(self.__dict__)

    def add_member(self, member):
        """
        Adds a LayoutDisplay to the database as part of this Layout
        """
        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
             INSERT INTO layout_displays (layout_id, display_id, grid_position, start_col, grid_width, start_row, grid_height)
                VALUES (%s, %s, %s, %s, %s, %s, %s)
        """
        cursor.execute(qry, (self.id, member.display_id, member.grid_position, member.start_col,
                             member.grid_width, member.start_row, member.grid_height))
        member.id = cursor.lastrowid
        self.members.append(member)

        cursor.close()
        conn.commit()
        conn.close()

    def dataset_ids(self):
        """
        Returns a list of the unique dataset IDs belonging to this layout
        """
        ids = set()

        for ds in self.members:
            ds.get_dataset_id()
            ids.add(ds.dataset_id)
        return list(ids)

    def display_ids(self):
        """
        Returns a list of the unique display IDs belonging to this layout
        """
        ids = set()

        for ds in self.members:
            if ds.display_id not in ids:
                ids.add(ds.display_id)

        return list(ids)

    def get_members(self, scope="all"):
        """
        Gets all members from the database and populates the 'members' attribute as a
        list of LayoutDisplay objects

        scope: "all" -> get all members, "single" -> get single-gene displays, "multi" -> get multi-gene displays

        """
        conn = Connection()
        cursor = conn.get_cursor()

        self.members = list()

        if scope not in ['all', 'single', 'multi']:
            raise Exception("ERROR: Invalid scope '{0}' passed to Layout.get_members()".format(scope))

        qry = """
            SELECT lm.id, lm.display_id, lm.grid_position, lm.start_col, lm.grid_width, lm.start_row, lm.grid_height
              FROM layout_displays lm
                   JOIN dataset_display d ON lm.display_id=d.id
                   JOIN dataset ds ON d.dataset_id=ds.id
             WHERE lm.layout_id = %s
                   AND ds.marked_for_removal = 0
        """

        cursor.execute(qry, (self.id,))

        for row in cursor:
            lm = LayoutDisplay(id=row[0], display_id=row[1], grid_position=row[2], start_col=row[3],
                                 grid_width=row[4], start_row=row[5], grid_height=row[6]
                )

            lm.get_dataset_id()
            lm.get_is_multigene()

            if lm.is_multigene:
                self.multigene_members.append(lm)
            else:
                self.singlegene_members.append(lm)

            if scope == "single" and lm.is_multigene:
                continue
            if scope == "multi" and not lm.is_multigene:
                continue

            self.members.append(lm)
        cursor.close()
        conn.close()

    def get_singlegene_members(self):
        return self.get_members(scope="single")

    def get_multigene_members(self):
        return self.get_members(scope="multi")

    def load(self):
        """
        If you only have the ID of layout and create an object from it, this method loads
        all the rest of the attributes, including layout members.
        """
        self.members = list()

        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              SELECT user_id, label, is_current, is_domain, share_id
                FROM layout
               WHERE id = %s
        """
        cursor.execute(qry, (self.id,))

        for row in cursor:
            (self.user_id, self.label, self.is_current, self.is_domain, self.share_id) = row

        self.get_members()

        cursor.close()
        conn.close()

    def remove(self):
        """
        Deletes the current layout from the database, along with its layout members.
        """
        self.remove_all_members()

        conn = Connection()
        cursor = conn.get_cursor()

        # First remove this as any user's default so it doesn't cause issues
        qry = "UPDATE guser SET layout_id = NULL WHERE layout_id = %s"
        cursor.execute(qry, (self.id,))

        qry = """
              DELETE FROM layout
              WHERE id = %s
        """
        cursor.execute(qry, (self.id,))

        cursor.close()
        conn.commit()
        conn.close()

    def remove_all_members(self):
        """
        Deletes all members for this layout from the database in one query
        """
        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              DELETE FROM layout_displays
              WHERE layout_id = %s
        """
        cursor.execute(qry, (self.id,))

        cursor.close()
        conn.commit()
        conn.close()

        self.members = []

    def remove_member_by_display_id(self, display_id):
        """Deletes a member for a given display ID from the database."""

        conn = Connection()
        cursor = conn.get_cursor()

        # limit removal to 1 member. TIe-breaker is highest grid_position
        qry = """
              DELETE FROM layout_displays
              WHERE display_id = %s
                AND layout_id = %s
                ORDER BY grid_position DESC LIMIT 1
        """
        cursor.execute(qry, (display_id, self.id))

        self.members = self.get_members()

        # ? Would be nice to delete empty rows, but then we have to figure out
        #  if this was a single or multigene display, then adjust the start_row for that type

        cursor.close()
        conn.commit()
        conn.close()

    def remove_members_by_dataset_id(self, dataset_id):
        """Deletes all members where the display ID belongs to a given dataset ID from the database."""

        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              DELETE FROM layout_displays
              WHERE display_id IN (SELECT
                                    id
                                    FROM dataset_display
                                    WHERE dataset_id = %s)
                AND layout_id = %s
        """

        cursor.execute(qry, (dataset_id, self.id))

        self.members = self.get_members()

        # ? Would be nice to dynamically adjust the start_row afterwards to remove gaps

        cursor.close()
        conn.commit()


    def save(self):
        """
        Will perform a save or an update depending on whether the ID attribute is
        defined.  If ID is already present, it's assumed to be an existing row which
        needs to be updated.
        """
        conn = Connection()
        cursor = conn.get_cursor()

        if self.id is None:
            layout_insert_qry = """
            INSERT INTO layout (user_id, label, is_current, is_domain, share_id)
            VALUES (%s, %s, %s, %s, %s)
            """
            cursor.execute(layout_insert_qry, (self.user_id, self.label, self.is_current,
                                               self.is_domain, self.share_id))
            self.id = cursor.lastrowid
        else:
            # ID already populated

            # Update layout properties
            sql = """
                  UPDATE layout
                     SET user_id = %s,
                         label = %s,
                         is_current = %s,
                         is_domain = %s,
                         share_id = %s
                   WHERE id = %s
            """
            cursor.execute(sql, (
                self.user_id, self.label, self.is_current,
                self.is_domain, self.share_id, self.id
            ))

            # TODO: delete existing members, add current ones

        cursor.close()
        conn.commit()
        conn.close()

    def save_change(self, attribute=None, value=None):
        """
        Update a layout attribute, both in the object and the relational database
        """
        if self.id is None:
            raise Exception("Error: no layout id. Cannot save change.")
        if attribute is None:
            raise Exception("Error: no attribute given. Cannot save change.")

        ## quick sanitization of attribute
        attribute = re.sub('[^a-zA-Z0-9_]', '_', attribute)
        setattr(self, attribute, value)

        conn = Connection()
        cursor = conn.get_cursor()

        save_sql = """
            UPDATE layout
            SET {0} = %s
            WHERE id = %s
        """.format(attribute)
        cursor.execute(save_sql, (str(value), self.id))

        conn.commit()
        cursor.close()
        conn.close()

    # TODO: Need a function to take a DatasetCollection and populate
    #  information on it within a layout

@dataclass
class LayoutCollection:
    # keep an index of folder IDs and their parent-most root IDs (tree walk needed here)
    root_folder_idx: dict = field(default_factory=dict, repr=False)

    # a simple index of folders and their parent IDs (this is a simple db query)
    folder_idx: dict = field(default_factory=dict, repr=False)

    layouts: List[Layout] = field(default_factory=list)

    # should dataset-populating methods called (adds overhead if you only need layout names)
    include_datasets: bool = True

    # In this class many of the methods need db connections, and this can be costly. Let's keep one
    #  open while it's used
    _cnx = Connection()

    def __post_init__(self):
        if len(self.folder_idx) == 0:
            self._populate_folder_index()
            self._populate_root_folder_index()

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        # Called when json modules attempts to serialize
        return self.__dict__

    def _get_root_folder_id(self, folder_id):
        """
        Recursive function to drive to the parent-most folder ID of
        any folder in the tree.
        Assumes self.folder_idx and self.root_folder_idx have been populated
        """

        # if the entry in folder_idx has a value, it's not a root node, so recurse.
        if self.folder_idx[folder_id]:
            return self._get_root_folder_id(self.folder_idx[folder_id])
        else:
            # else this is a root node
            return folder_id

    def _populate_folder_index(self):
        """
        Populates the self.folder_idx attribute with a dictionary where the index is
        folder ID and value is that folder's parent ID.
        """
        self.folder_idx = dict()

        qry = "SELECT id, parent_id FROM folder"
        cursor = self._cnx.get_cursor()
        cursor.execute(qry)

        for row in cursor:
            self.folder_idx[row[0]] = row[1]

        cursor.close()

    def _populate_root_folder_index(self):
        """
        Populates the self.root_folder_idx attribute with a dictionary where the index is
        folder ID and value is root folder ID the ones in gear.ini[folders]
        """
        if len(self.folder_idx):
            return False

        qry = "SELECT id, parent_id FROM folder"
        cursor = self._cnx.get_cursor()
        cursor.execute(qry)

        for row in cursor:
            # if the parent_id is empty, we don't need to store it because it's a top-level node
            if row[1]:
                self.root_folder_idx[row[0]] = self._get_root_folder_id(row[1])

        cursor.close()

    def get_by_share_id(self, share_id=None):
        """
        Gets the layout from the passed share_id, if any.
        TODO: f.parent_id, f.label can be removed from all these get_by* methods
        """
        if not share_id:
            return self.layouts

        cursor = self._cnx.get_cursor()

        qry = """
              SELECT l.id, l.label, l.is_current, l.user_id, l.share_id, l.is_domain, l.is_public
                FROM layout l
               WHERE l.share_id = %s
        """
        cursor.execute(qry, (share_id,))

        for row in cursor:
            layout = Layout(
                id=row[0],
                label=row[1],
                is_current=row[2],
                user_id=row[3],
                share_id=row[4],
                is_domain=row[5],
                is_public=row[6]
            )


            if self.include_datasets == True:
                layout.get_members()
                layout.dataset_count = len(layout.members)  # Excludes datasets marked for removal

            self.layouts.append(layout)

        cursor.close()
        return self.layouts

    def get_by_user(self, user=None):
        """
        Gets all the layouts owned by a user
        """
        if not isinstance(user, User):
            raise Exception("LayoutCollection.get_by_user() requires an instance of User to be passed.")

        cursor = self._cnx.get_cursor()

        qry = """
              SELECT l.id, l.label, l.is_current, l.user_id, l.share_id, l.is_domain, l.is_public
                FROM layout l
               WHERE l.user_id = %s
        """
        cursor.execute(qry, (user.id,))

        for row in cursor:
            layout = Layout(
                id=row[0],
                label=row[1],
                is_current=row[2],
                user_id=row[3],
                share_id=row[4],
                is_domain=row[5],
                is_public=row[6]
            )

            if self.include_datasets == True:
                layout.get_members()
                layout.dataset_count = len(layout.members)  # Excludes datasets marked for removal
                print(layout.members, file=sys.stderr)

            self.layouts.append(layout)

        cursor.close()
        return self.layouts

    def get_by_users_groups(self, user=None, append=True):
        """
        Queries the DB to get all the groups of which the passed user is a member, then
        gets all layouts in those groups.
        If the append argument is set to False, the class' layouts attribute will be
        cleared before these are aded.
        """
        if not isinstance(user, User):
            raise Exception("LayoutCollection.get_by_users_groups() requires an instance of User to be passed.")

        if append == False:
            self.layouts = list()

        cursor = self._cnx.get_cursor()

        qry = """
              SELECT l.id, l.label, l.is_current, l.user_id, l.share_id, l.is_domain, l.is_public
                FROM ggroup g
                     JOIN user_group_membership ugm ON ugm.group_id=g.id
                     JOIN guser u ON u.id=ugm.user_id
                     JOIN layout_group_membership lgm ON lgm.group_id=g.id
                     JOIN layout l ON lgm.layout_id=l.id
               WHERE u.id = %s
        """
        cursor.execute(qry, (user.id,))

        for row in cursor:
            layout = Layout(
                id=row[0],
                label=row[1],
                is_current=row[2],
                user_id=row[3],
                share_id=row[4],
                is_domain=row[5],
                is_public=row[6],
            )

            if self.include_datasets == True:
                layout.get_members()
                layout.dataset_count = len(layout.members)  # Excludes datasets marked for removal

            self.layouts.append(layout)

        cursor.close()
        return self.layouts

    def get_domains(self):
        """
        Queries the DB to get all the site domain layouts.
        """
        cursor = self._cnx.get_cursor()

        qry = """
              SELECT l.id, l.label, l.is_current, l.user_id, l.share_id, l.is_domain, l.is_public
                FROM layout l
               WHERE l.is_domain = 1
        """

        cursor.execute(qry)

        for row in cursor:
            layout = Layout(
                id=row[0],
                label=row[1],
                is_current=row[2],
                user_id=row[3],
                share_id=row[4],
                is_domain=row[5],
                is_public=row[6]
            )

            if self.include_datasets == True:
                layout.get_members()
                layout.dataset_count = len(layout.members)  # Excludes datasets marked for removal

            self.layouts.append(layout)

        cursor.close()
        return self.layouts

    def get_public(self):
        """
        Queries the DB to get all the public layouts.
        """
        cursor = self._cnx.get_cursor()

        qry = """
              SELECT l.id, l.label, l.is_current, l.user_id, l.share_id, l.is_domain, l.is_public
                FROM layout l
               WHERE l.is_public = 1
        """

        cursor.execute(qry)

        for row in cursor:
            layout = Layout(
                id=row[0],
                label=row[1],
                is_current=row[2],
                user_id=row[3],
                share_id=row[4],
                is_domain=row[5],
                is_public=row[6]
            )

            if self.include_datasets == True:
                layout.get_members()
                layout.dataset_count = len(layout.members)  # Excludes datasets marked for removal

            self.layouts.append(layout)

        cursor.close()
        return self.layouts

@dataclass
class Folder:
    id: int = None
    parent_id: int = None
    label: str = None

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        return self.__dict__

@dataclass
class FolderCollection:
    folders: List[Folder] = field(default_factory=list)

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        # Called when json modules attempts to serialize
        return self.__dict__

    def get_by_folder_ids(self, ids=None):
        """
        Populates the self.folders attribute with the passed list of
        folder IDs. Resets self.folders on each execution.
        """
        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)
        self.folders = list()

        ## Sanitize the IDs
        cleaned = [ str(x) for x in ids if isinstance(x, int) ]

        qry = """
              SELECT id, parent_id, label
                FROM folder
               WHERE id in ({0})
        """.format(",".join(cleaned))

        cursor.execute(qry)

        for row in cursor:
            folder = Folder(row)
            self.folders.append(folder)

        cursor.close()
        conn.close()
        return self.folders

    def get_root_folders(self, folder_type=None):
        """
        Returns a list of Folder elements for the root folders.
        These should all map to entries within gear.ini -> [folders]
        folder_type should be either 'cart' or 'profile'
        """
        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = "SELECT id, label FROM folder WHERE parent_id IS NULL";
        cursor.execute(qry)

        # Create an index of all the root folders.  Still need to filter/map
        # them based on their type
        rfs = dict()
        for row in cursor:
            rfs[row['id']] = row['label']

        for scope in ['domain', 'user', 'group', 'shared', 'public']:
            config_key = "{0}_{1}_master_id".format(folder_type, scope)
            folder_id = int(this.servercfg['folders'][config_key])

            if folder_id in rfs:
                folder = Folder(id=folder_id,
                                parent_id=None,
                                label=rfs[folder_id])
                self.folders.append(folder)

        cursor.close()
        conn.close()
        return self.folders

    def get_tree_by_folder_ids(self, ids=None, folder_type=None):
        """
        Similar to get_by_folder_ids() but this gets the entire tree for any folder IDs
        passed. Higher query cost than just running get_by_folder_ids().
        Also always returns the root folders, even if they're empty.
        This link looked like a good solution to handle within the database directly but
        couldn't get it to work across MySQL/MariaDB and versions.  Too delicate:
        https://stackoverflow.com/a/60019201/1368079
        """
        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)
        self.get_root_folders(folder_type=folder_type)
        all_ids = list()

        # Remember the root folders so we don't duplicate them
        for folder in self.folders:
            all_ids.append(folder.id)

        new_ids_found = ids

        while len(new_ids_found) > 0:
            cleaned = [ str(x) for x in new_ids_found if isinstance(x, int) ]

            qry = """
                  SELECT id, parent_id, label
                    FROM folder
                   WHERE id in ({0})
            """.format(",".join(cleaned))

            cursor.execute(qry)
            new_ids_found = list()

            for row in cursor:
                folder = Folder(id=row['id'], parent_id=row['parent_id'], label=row['label'])

                if folder.id not in all_ids:
                    all_ids.append(folder.id)
                    self.folders.append(folder)

                if folder.parent_id and folder.parent_id not in all_ids:
                    new_ids_found.append(folder.parent_id)

        cursor.close()
        conn.close()
        return self.folders

@dataclass
class DatasetLink:
    id: Optional[int] = None
    dataset_id: Optional[str] = None
    resource: Optional[str] = None
    label: Optional[str] = None
    url: Optional[str] = None

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        return self.__dict__

@dataclass
class DatasetDisplay:
    id: Optional[int] = None
    dataset_id: Optional[str] = None
    user_id: Optional[int] = None
    label: Optional[str] = None
    plot_type: Optional[str] = None
    plotly_config: Optional[str] = None

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        return self.__dict__

    def save(self):
        """
        Currently assumes the display already exists and what you're saving here is a change
        to it.
        """
        if self.id is None:
            # TODO: implement new saves
            raise Exception("Can only save changes to existing dataset displays currently")
        else:
            # gene_cart_member entries are deleted by foreign key cascade
            conn = Connection()
            cursor = conn.get_cursor()

            sql = """
                  UPDATE dataset_display
                     SET dataset_id = %s,
                         user_id = %s,
                         label = %s,
                         plot_type = %s,
                         plotly_config = %s
                   WHERE id = %s
            """
            cursor.execute(sql, (
                self.dataset_id, self.user_id, self.label,
                self.plot_type, self.plotly_config, self.id
            ))

            conn.commit()
            cursor.close()
            conn.close()

    def remove(self):
        """
        Deletes the current display from the database.
        """
        conn = Connection()
        cursor = conn.get_cursor()

        # first remove any layout_display entries
        qry = """
                DELETE FROM layout_displays
                WHERE display_id = %s
        """
        cursor.execute(qry, (self.id,))

        # Then remove from anly dataset_preference entries
        qry = """
                DELETE FROM dataset_preference
                WHERE display_id = %s
        """
        cursor.execute(qry, (self.id,))

        # Then remove the display itself
        qry = """
                DELETE FROM dataset_display
                WHERE id = %s
        """
        cursor.execute(qry, (self.id,))

        cursor.close()
        conn.commit()
        conn.close()

@dataclass
class Dataset:
    id: str
    owner_id: Optional[int] = None
    title: Optional[str] = None
    organism_id: Optional[int] = None
    pubmed_id: Optional[str] = None
    geo_id: Optional[str] = None
    is_public: Optional[int] = None
    is_downloadable: Optional[int] = None
    ldesc: Optional[str] = None
    date_added: Optional[datetime.datetime] = None
    dtype: Optional[str] = None
    schematic_image: Optional[str] = None
    share_id: Optional[str] = None
    math_default: Optional[str] = None
    marked_for_removal: Optional[int] = None
    load_status: Optional[str] = None
    has_h5ad: Optional[int] = None
    platform_id: Optional[str] = None
    instrument_model: Optional[str] = None
    library_selection: Optional[str] = None
    library_source: Optional[str] = None
    library_strategy: Optional[str] = None
    contact_email: Optional[str] = None
    contact_institute: Optional[str] = None
    contact_name: Optional[str] = None
    annotation_source: Optional[str] = None
    plot_default: Optional[str] = None
    annotation_release: Optional[int] = None
    # derived, here for convenience
    gene_count: Optional[int] = None
    obs_count: Optional[int] = None
    has_tarball: int = 0
    displays: List[DatasetDisplay] = field(default_factory=list)
    tags: List[str] = field(default_factory=list)
    layouts: List[Layout] = field(default_factory=list)
    links: List[DatasetLink] = field(default_factory=list)
    user_name: Optional[str] = None
    is_permalink: Optional[int] = None

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        return self.__dict__

    def get_displays(self):
        """
        Populates the dataset displays attribute, a list of DatasetDisplay objects
        related to this dataset.
        """
        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              SELECT id, user_id, label, plot_type, plotly_config
                FROM dataset_display
               WHERE dataset_id = %s
        """
        cursor.execute(qry, (self.id,))
        self.displays.clear()

        for (display_id, user_id, label, plot_type, plotly_config) in cursor:
            display = DatasetDisplay(
                id=display_id,
                dataset_id=self.id,
                user_id=user_id,
                label=label,
                plot_type=plot_type,
                plotly_config=plotly_config
            )

            self.displays.append(display)

        cursor.close()
        conn.close()

        return self.displays

    def get_owner_display(self, is_multigene=False):
        """
        Returns the display object for the owner of the dataset.
        """
        if self.owner_id is None:
            return None

        qry = """
              SELECT dd.id, dd.label, dd.plot_type, dd.plotly_config
                FROM dataset d
                    JOIN dataset_display dd ON d.id = dd.dataset_id
                    JOIN dataset_preference dp ON d.id=dp.dataset_id
                WHERE dp.user_id = d.owner_id
                AND dd.user_id = d.owner_id
                AND dp.is_multigene = %s
                AND d.id = %s;
        """
        conn = Connection()
        cursor = conn.get_cursor()

        cursor.execute(qry, (is_multigene, self.id))

        for row in cursor:
            display = DatasetDisplay(
                id=row[0],
                dataset_id=self.id,
                user_id=self.owner_id,
                label=row[1],
                plot_type=row[2],
                plotly_config=row[3]
            )

            cursor.close()
            conn.close()
            return display

        return None

    def get_file_path(self, session_id=None):
        """
        The file path of a dataset can depend on the dataset type as well as whether we
        are looking at the primary file or one generated by a user session during an
        analysis.

        If a session_id is passed, it's assumed that we'll look for a session-specific copy.
        If one isn't found, the original source file is returned.

        This returns where the path SHOULD be, it doesn't check that it's actually there. This
        allows for it to be used also for any process which wants to know where to write it.
        """

        if self.dtype == "spatial":
            if session_id is None:
                zarr_file_path = "{0}/../www/datasets/spatial/{1}.zarr".format(
                    os.path.dirname(os.path.abspath(__file__)), self.id)
            else:
                zarr_file_path = "{0}/{1}/spatial/{2}.zarr".format(this.analysis_base_dir, session_id, self.id)
            return zarr_file_path

        if session_id is None:
            h5ad_file_path = "{0}/../www/datasets/{1}.h5ad".format(
                os.path.dirname(os.path.abspath(__file__)), self.id)
        else:
            h5ad_file_path = "{0}/{1}/{2}.h5ad".format(this.analysis_base_dir, session_id, self.id)

        return h5ad_file_path

    def get_tarball_path(self):
        """
        This returns where the path of where a dataset's tarball SHOULD be, it doesn't check that
        it's actually there. This allows for it to be used also for any process which wants to
        know where to write it.
        """

        tarball_file_path = "{0}/../www/datasets/{1}.tar.gz".format(
            os.path.dirname(os.path.abspath(__file__)), self.id)

        return tarball_file_path

    def get_layouts(self, user=None, include_public=False):
        """
        Populates the dataset layouts attribute, a list of Layout objects in which
        this dataset can be found (only those which the user has rights to see.)

        First checks to see if self.layouts is empty.  If already populated, it is
        just returned.
        """
        if len(self.layouts) < 1:
            conn = Connection()
            cursor = conn.get_cursor()

            if user:
                qry = """
                      SELECT DISTINCT l.id, l.user_id, l.is_domain, l.label, l.is_current, l.share_id
                        FROM layout l
                             JOIN layout_displays lm ON lm.layout_id=l.id
                             JOIN dataset_display dd ON dd.id =lm.display_id
                       WHERE dd.dataset_id = %s
                             AND (l.is_domain = 1 OR l.user_id = %s)
                    ORDER BY l.label
                """
                cursor.execute(qry, (self.id, user.id))
            else:
                qry = """
                      SELECT DISTINCT l.id, l.user_id, l.is_domain, l.label, l.is_current, l.share_id
                        FROM layout l
                             JOIN layout_displays lm ON lm.layout_id=l.id
                             JOIN dataset_display dd ON dd.id =lm.display_id
                       WHERE dd.dataset_id = %s
                             AND l.is_domain = 1
                    ORDER BY l.label
                """
                cursor.execute(qry, (self.id,))

            for row in cursor:
                l = Layout(id=row[0], user_id=row[1], is_domain=row[2], label=row[3],
                           is_current=row[4], share_id=row[5])
                self.layouts.append(l)

            if include_public:
                qry = """
                      SELECT DISTINCT l.id, l.user_id, l.is_domain, l.label, l.is_current, l.share_id
                        FROM layout l
                             JOIN layout_displays lm ON lm.layout_id=l.id
                             JOIN dataset_display dd ON dd.id =lm.display_id
                       WHERE dd.dataset_id = %s
                             AND l.is_public = 1
                             AND l.is_domain = 0
                    ORDER BY l.label
                """
                cursor.execute(qry, (self.id,))

                for row in cursor:
                    l = Layout(id=row[0], user_id=row[1], is_domain=row[2], label=row[3],
                               is_current=row[4], share_id=row[5])
                    self.layouts.append(l)

            # deduplicate based on a layout's ID (in case user and public layouts overlap)
            # (normal set->list conversion doesn't work for objects)
            self.layouts = list({l.id: l for l in self.layouts}.values())

            cursor.close()

        return self.layouts

    def get_links(self):
        """
        Populates the dataset links attribute, a list of DatasetLink objects
        associated with this dataset.

        First checks to see if self.links is empty.  If already populated, it is
        just returned.
        """
        if len(self.links) < 1:
            conn = Connection()
            cursor = conn.get_cursor()

            qry = "SELECT id, resource, label, url FROM dataset_link WHERE dataset_id = %s"
            cursor.execute(qry, (self.id,))

            for (id, resource, label, url) in cursor:
                dsl = DatasetLink(dataset_id=self.id, resource=resource, label=label, url=url)
                self.links.append(dsl)

            cursor.close()

        return self.links

    def get_shape(self, session_id=None, tuple_only=False):
        """
        Queries the dataset's source expression matrix in order to get its shape.

        This updates the gene_count and obs_count attributes, returns a string
        like '10000x20000' in format genes x obs
        """
        if self.dtype is None:
            raise Exception("Error: can't get shape for a dataset of unknown type")
        elif self.has_h5ad:
            ## File is under datasets/${id}.h5ad
            h5ad_file_path = self.get_file_path(session_id=session_id)

            from shadows import AnnDataShadow
            adata = AnnDataShadow(h5ad_file_path)

            (n_obs, n_vars) = adata.shape

            adata.close()

            if tuple_only:
                return (n_obs, n_vars)

            self.gene_count = n_vars
            self.obs_count = n_obs

            return "{0}x{1}".format(self.gene_count, self.obs_count)

    def to_json(self):
        return str(self)

    def remove(self):
        """
        Removes a dataset and its dependencies from the database
        """
        raise Exception("Support not yet added to remove dataset via the API")

    def save_change(self, attribute=None, value=None):
        """
        Update a dataset attribute, both in the object and the relational database
        """
        if self.id is None:
            raise Exception("Error: no dataset id. Cannot save change.")
        if attribute is None:
            raise Exception("Error: no attribute given. Cannot save change.")

        ## quick sanitization of attribute
        attribute = re.sub('[^a-zA-Z0-9_]', '_', attribute)
        setattr(self, attribute, value)

        conn = Connection()
        cursor = conn.get_cursor()

        save_sql = """
            UPDATE dataset
            SET {0} = %s
            WHERE id = %s
        """.format(attribute)
        cursor.execute(save_sql, (str(value), self.id))

        conn.commit()
        cursor.close()
        conn.close()


@dataclass
class DatasetCollection:
    datasets: List[Dataset] = field(default_factory=list)

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        # Called when json modules attempts to serialize
        return self.__dict__

    def filter_by_types(self, types=None):
        """
        Filters the collection by one or more dataset types (a list).
        """
        datasets_to_keep = []

        for dataset in self.datasets:
            if dataset.dtype in types:
                datasets_to_keep.append(dataset)

        self.datasets = datasets_to_keep

    def get_by_dataset_ids(self, ids=None, get_links=False):
        conn = Connection()
        cursor = conn.get_cursor()

        # I tried but couldn't make this work in a single query with IN () instead.
        qry = """
           SELECT d.id, d.title, o.label, d.pubmed_id, d.geo_id, d.is_public, d.is_downloadable, d.ldesc,
                  d.dtype, u.id, u.user_name, d.schematic_image, d.share_id,
                  d.math_default, d.marked_for_removal, d.date_added, d.load_status,
                  IFNULL(GROUP_CONCAT(t.label), 'NULL') as tags,
                  d.annotation_source, d.annotation_release, d.organism_id
             FROM dataset d
                  JOIN organism o ON d.organism_id=o.id
                  JOIN guser u ON d.owner_id=u.id
                  LEFT JOIN dataset_tag dt ON dt.dataset_id = IFNULL(d.id, 'NULL')
                  LEFT JOIN tag t ON t.id = IFNULL(dt.tag_id, 'NULL')
            WHERE d.id = %s
              AND d.marked_for_removal != 1
        GROUP BY d.id, d.title, o.label, d.pubmed_id, d.geo_id, d.is_public, d.is_downloadable, d.ldesc,
                   d.dtype, u.id, u.user_name, d.schematic_image, d.share_id,
                   d.math_default, d.marked_for_removal, d.date_added, d.load_status,
                   d.annotation_source, d.annotation_release, d.organism_id
        """
        for id in ids:
            cursor.execute(qry, (id,))

            for row in cursor:

                # valid pubmed IDs are numeric
                m = re.match(r"\d+", str(row[3]))
                if m:
                    pubmed_id = row[3]
                else:
                    pubmed_id = None

                date_added = row[15].isoformat()

                dataset = Dataset(id=row[0],
                                title=row[1],
                                organism_id=row[20],
                                pubmed_id=pubmed_id,
                                geo_id=row[4],
                                is_public=row[5],
                                is_downloadable=row[6],
                                ldesc=row[7],
                                dtype=row[8],
                                owner_id=row[9],
                                schematic_image=row[11],
                                share_id=row[12],
                                math_default=row[13],
                                date_added=date_added,
                                load_status=row[16],
                                annotation_source=row[18],
                                annotation_release=row[19]
                )

                # Add supplemental attributes this method created previously
                dataset.organism = row[2]
                dataset.tags = row[17].split(',')
                dataset.access = 'access_level'
                dataset.user_name = row[10]

                if dataset.dtype == "spatial":
                    #NotImplementedError("Spatial datasets don't have permanent tarball location yet")
                    dataset.has_tarball = 0
                    dataset.has_h5ad = 0
                else:
                    if os.path.exists(dataset.get_tarball_path()):
                        dataset.has_tarball = 1
                    else:
                        dataset.has_tarball = 0

                    if os.path.exists(dataset.get_file_path()):
                        dataset.has_h5ad = 1
                    else:
                        dataset.has_h5ad = 0

                #  TODO: These all need to be tracked through the code and removed
                dataset.dataset_id = dataset.id
                dataset.user_id = dataset.owner_id
                dataset.math_format = dataset.math_default
                # plot_default needs to be removed everywhere
                # load_status needs to be removed everywhere
                # schematic_image ?

                if get_links:
                    dataset.get_links()

                self.datasets.append(dataset)

        cursor.close()
        conn.close()
        return self.datasets

    def get_public(self, has_h5ad=None, n=None, order_by=None, types=None):
        """
        Populates the DatasetCollection's datasets list with those datasets which are marked
        as public (with optional filtering arguments.)

        Does NOT clear the existing internal list first, so a few methods could be called to
        make custom collections.

        order_by: ['title' (default), 'date_added']

        types: ['single-cell-rnaseq', 'microarray']

        If ordering by date_added, the newest will be first.
        """

        if order_by is None:
            order_by = 'title'

        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = """
              SELECT id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, is_downloadable, ldesc,
                     date_added, dtype, schematic_image, share_id, math_default, marked_for_removal,
                     load_status, has_h5ad
                FROM dataset d
               WHERE d.is_public = 1
                 AND d.marked_for_removal = 0
        """
        qry_args = []

        if has_h5ad is not None:
            qry += " AND d.has_h5ad = %s"
            qry_args.append(has_h5ad)

        if order_by == 'title':
            qry += " ORDER BY d.title"
        elif order_by == 'date_added':
            qry += " ORDER BY d.date_added DESC"

        if n is not None:
            qry += " LIMIT {0}".format(n)

        cursor.execute(qry, qry_args)
        for r in cursor:
            dataset = Dataset(id=r['id'], owner_id=r['owner_id'], title=r['title'],
                              organism_id=r['organism_id'], pubmed_id=r['pubmed_id'],
                              geo_id=r['geo_id'], is_public=r['is_public'], is_downloadable=r['is_downloadable'], ldesc=r['ldesc'],
                              date_added=r['date_added'], dtype=r['dtype'],
                              schematic_image=r['schematic_image'], share_id=r['share_id'],
                              math_default=r['math_default'], marked_for_removal=r['marked_for_removal'],
                              load_status=r['load_status'], has_h5ad=r['has_h5ad'])
            self.datasets.append(dataset)

        if types is not None:
            self.filter_by_types(types=types)

        cursor.close()
        conn.close()

    def get_owned_by_user(self, has_h5ad=None, user=None, types=None):
        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = """
              SELECT d.id, d.owner_id, d.title, organism_id, pubmed_id, geo_id, is_public, is_downloadable, ldesc,
                     date_added, dtype, schematic_image, share_id, math_default, marked_for_removal,
                     load_status, has_h5ad
                FROM dataset d
               WHERE d.owner_id = %s
                 AND d.marked_for_removal = 0
        """
        qry_args = [user.id]

        if has_h5ad is not None:
            qry += " AND d.has_h5ad = %s"
            qry_args.append(has_h5ad)

        qry += " ORDER BY d.title"

        cursor.execute(qry, qry_args)
        for r in cursor:
            dataset = Dataset(id=r['id'], owner_id=r['owner_id'], title=r['title'],
                              organism_id=r['organism_id'], pubmed_id=r['pubmed_id'],
                              geo_id=r['geo_id'], is_public=r['is_public'], is_downloadable=r['is_downloadable'], ldesc=r['ldesc'],
                              date_added=r['date_added'], dtype=r['dtype'],
                              schematic_image=r['schematic_image'], share_id=r['share_id'],
                              math_default=r['math_default'], marked_for_removal=r['marked_for_removal'],
                              load_status=r['load_status'], has_h5ad=r['has_h5ad'])
            self.datasets.append(dataset)

        if types is not None:
            self.filter_by_types(types=types)

        cursor.close()
        conn.close()


    def get_shared_with_user(self, has_h5ad=None, user=None, types=None):
        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = """
              SELECT d.id, d.owner_id, d.title, organism_id, pubmed_id, geo_id, is_public, is_downloadable, ldesc,
                     date_added, dtype, schematic_image, share_id, math_default, marked_for_removal,
                     load_status, has_h5ad
                FROM dataset d
                     JOIN dataset_shares ds ON d.id=ds.dataset_id
               WHERE ds.user_id = %s
                 AND ds.is_allowed = 1
                 AND d.marked_for_removal = 0
        """
        qry_args = [user.id]

        if has_h5ad is not None:
            qry += " AND d.has_h5ad = %s"
            qry_args.append(has_h5ad)

        qry += " ORDER BY d.title"

        cursor.execute(qry, qry_args)
        for r in cursor:
            dataset = Dataset(id=r['id'], owner_id=r['owner_id'], title=r['title'],
                              organism_id=r['organism_id'], pubmed_id=r['pubmed_id'],
                              geo_id=r['geo_id'], is_public=r['is_public'], is_downloadable=r['is_downloadable'], ldesc=r['ldesc'],
                              date_added=r['date_added'], dtype=r['dtype'],
                              schematic_image=r['schematic_image'], share_id=r['share_id'],
                              math_default=r['math_default'], marked_for_removal=r['marked_for_removal'],
                              load_status=r['load_status'], has_h5ad=r['has_h5ad'])
            self.datasets.append(dataset)

        if types is not None:
            self.filter_by_types(types=types)

        cursor.close()
        conn.close()

    def to_json(self):
        return str(self)

@dataclass
class Gene:
    id: Optional[int] = None
    ensembl_id: Optional[str] = None
    ensembl_version: Optional[str] = None
    ensembl_release: Optional[int] = None
    genbank_acc: Optional[str] = None
    organism_id: Optional[int] = None
    molecule: Optional[str] = None
    start: Optional[int] = None
    stop: Optional[int] = None
    gene_symbol: Optional[str] = None
    product: Optional[str] = None
    biotype: Optional[str] = None

    # derived, not in the relational DB
    go_terms: List[dict] = field(default_factory=list)
    dbxrefs: List[dict] = field(default_factory=list)
    aliases: List[dict] = field(default_factory=list)

    def __repr__(self):
        return json.dumps(self.__dict__)

    def load_aliases(self):
        self.aliases = list()

        qry = "SELECT label FROM gene_symbol WHERE gene_id = %s AND is_primary = 0"
        conn = Connection()
        cursor = conn.get_cursor()
        cursor.execute(qry, (self.id, ))

        for (alias,) in cursor:
            self.aliases.append({'label': alias})

        cursor.close()
        conn.close()

    def load_dbxrefs(self):
        self.dbxrefs = list()

        qry = "SELECT dbxref FROM gene_dbxref WHERE gene_id = %s"
        conn = Connection()
        cursor = conn.get_cursor()
        cursor.execute(qry, (self.id, ))

        for (dbxref,) in cursor:
            dbxref_parts = dbxref.split(':')
            source = dbxref_parts[0]
            identifier = dbxref_parts[1]
            self.dbxrefs.append({'source': source, 'identifier': identifier, 'url': None})

        cursor.close()
        conn.close()

    def load_dbxref_links(self):
        """
        Loads the database references for the current gene and also forms URL strings for
        each where we have a mapping from the source to a URL.  If there are no current
        dbxrefs it pulls them first.  If some are present, it trusts that is the list
        to use and doesn't query the database for them.

        Each link element has these attributes:

        - source: The label which will show
        - identifier: What will be used in URL building, usually the gene symbol
        - url: Well, the URL
        - title: What shows on hover, often an expansion of the source (optional)
        """
        if len(self.dbxrefs) == 0:
            self.load_dbxrefs()

        # Do the Ensembl one
        self.dbxrefs.append({'source': 'ENSEMBL', 'identifier': self.ensembl_id,
                             'url': "https://www.ensembl.org/id/{0}".format(self.ensembl_id)})

        # Add a PubMed search term link
        self.dbxrefs.append({'source': 'PubMed', 'identifier': self.gene_symbol,
                             'url': "http://www.ncbi.nlm.nih.gov/pubmed/?term={0}".format(self.gene_symbol)})

        # Does the user want Homologene links?
        if 'HomoloGene' in this.links_out and this.links_out['HomoloGene'] is True:
            hg_url = 'https://www.ncbi.nlm.nih.gov/homologene/?term=' + self.gene_symbol.upper()
            self.dbxrefs.append({'source': 'HomoloGene', 'identifier': self.gene_symbol, 'url': hg_url})

        if this.domain_label == 'gear':
            # Add one for the SHIELD
            self.dbxrefs.append({'source': 'SHIELD', 'identifier': self.gene_symbol, 'url': None})

            # Hack to add DVD links until they have SSL enabled and we can use their API instead
            if 'DVD' in this.links_out and this.links_out['DVD'] is True:
                dvd_genes = ["ACTG1","ADCY1","ADGRV1","GPR98","AIFM1","ALMS1","ATP2B2","ATP6V1B1","BDP1","BSND","C10orf2","CABP2","CACNA1D","CCDC50","CD164","CDC14A","CDH23","CEACAM16","CIB2","CISD2","CLDN14","CLIC5","CLPP","CLRN1","COCH","COL11A1","COL11A2","COL2A1","COL4A3","COL4A4","COL4A5","COL4A6","COL9A1","COL9A2","CRYM","DCDC2","DFNA5","DFNB31 (WHRN)","DFNB59 (PJVK)","DIABLO","DIAPH1","DIAPH3","DSPP","EDN3","EDNRB","ELMOD3","EPS8","EPS8L2","ESPN","ESRRB","EYA1","EYA4","FAM65B","FGF3","FGFR1","FGFR2","FOXI1","GATA3","GIPC3","GJB2","GJB3","GJB6","GPSM2","GRHL2","GRXCR1","GRXCR2","HARS2","HGF","HOMER2","HSD17B4","ILDR1","KARS","KCNE1","KCNJ10","KCNQ1","KCNQ4","KITLG","LARS2","LHFPL5","LOXHD1","LOXL3","LRTOMT","MARVELD2","MCM2","MET","MIR96","MITF","MSRB3","MT-RNR1","MT-TL1","MT-TS1","MYH14","MYH9","MYO15A","MYO3A","MYO6","MYO7A","NARS2","NLRP3","OPA1","OSBPL2","OTOA","OTOF","OTOG","OTOGL","P2RX2","PAX3","PCDH15","PDZD7","PEX1","PEX6","PNPT1","POLR1C","POLR1D","POU3F4","POU4F3","PRPS1","PTPRQ","RDX","ROR1","S1PR2","SERPINB6","SIX1","SIX5","SLC17A8","SLC22A4","SLC26A4","SLC26A5","SLITRK6","SMPX","SNAI2","SOX10","STRC","SYNE4","TBC1D24","TBX1","TCOF1","TECTA","TECTB","TIMM8A","TJP2","TMC1","TMEM132E","TMIE","TMPRSS3","TNC","TPRN","TRIOBP","TSPEAR","USH1C","USH1G","USH2A","WFS1"]
                if self.gene_symbol.upper() in dvd_genes:
                    dvd_url = 'http://deafnessvariationdatabase.org/gene/' + self.gene_symbol.upper()
                    self.dbxrefs.append({'source': 'DVD', 'identifier': self.gene_symbol, 'url': dvd_url, 'title': 'Deafness Variation Database'})

        elif this.domain_label == 'nemo':
            pass

        for dbxref in self.dbxrefs:
            source = dbxref['source']

            if source == 'CCDS':
                dbxref['url'] = "https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&GO=MainBrowse&DATA={0}".format(dbxref['identifier'])
            elif source == 'SHIELD':
                dbxref['url'] = "https://shield.hms.harvard.edu/viewgene.html?gene={0}".format(dbxref['identifier'])
            elif source == 'UCSC':
                if self.organism_id == 1:
                    dbxref['url'] = ("http://genome.ucsc.edu/cgi-bin/hgTracks?org=mouse&db=mm10&"
                                     "singleSearch=knownCanonical&position={0}".format(dbxref['identifier']))
            elif source == 'UniParc':
                dbxref['url'] = "https://www.uniprot.org/uniparc/{0}?sort=score".format(dbxref['identifier'])
            elif source == 'Uniprot/SWISSPROT':
                dbxref['url'] = "https://www.uniprot.org/uniprot/{0}".format(dbxref['identifier'])

                # Bonus, all swissprot entries are also indexed in Phosphosite
                phospho_url = "https://www.phosphosite.org/uniprotAccAction?id={0}".format(dbxref['identifier'])
                self.dbxrefs.append({'source': 'PhosphoSite', 'identifier': self.gene_symbol, 'url': phospho_url})

    def load_go_terms(self):
        self.go_terms = list()

        qry = """
              SELECT go.go_id, go.name
               FROM gene_go_link ggl
                    JOIN go ON ggl.go_id = go.go_id
               WHERE ggl.gene_id = %s
        """
        conn = Connection()
        cursor = conn.get_cursor()
        cursor.execute(qry, (self.id, ))

        for (go_id, go_name) in cursor:
            self.go_terms.append({'go_id':go_id, 'name':go_name})

        cursor.close()
        conn.close()



@dataclass
class GeneCollection:
    genes: List[Gene] = field(default_factory=list)

    def add_gene(self, gene):
        self.genes.append(gene)

    def get_by_gene_symbol(self, gene_symbol=None, exact=None, organism_id=None):
        """
        Searches the database by gene symbol and populates the GeneCollection.genes attribute.

        The 'gene_symbol' argument can be a space-separated list
        """
        conn = Connection()
        cursor = conn.get_cursor()

        # Keeps track of those we've found already so we don't show duplicates.
        #  The count is how many of that gene symbol were in that release.
        # h{gene_symbol}{organism_id} = {'ensembl_release': 93, count: 3}
        found = defaultdict(dict)

        gene_symbols = gene_symbol.split(' ')

        org_addendum = ""
        if organism_id:
            org_addendum = "AND organism_id = {} ".format(organism_id)

        for gene_symbol in gene_symbols:
            if len(gene_symbol) == 0:
                continue

            # this comes from javascript as a string true/false
            if exact in [True, 'true', 1, "1"]:
                qry = """
                  SELECT id, ensembl_id, ensembl_version, ensembl_release, genbank_acc, organism_id,
                         molecule, start, stop, gene_symbol, product, biotype
                    FROM gene
                   WHERE gene_symbol = %s {}
                  ORDER BY gene_symbol, organism_id, ensembl_release DESC
                """.format(org_addendum)
                cursor.execute(qry, (gene_symbol,))
            else:
                qry = """
                  SELECT id, ensembl_id, ensembl_version, ensembl_release, genbank_acc, organism_id,
                         molecule, start, stop, gene_symbol, product, biotype
                    FROM gene
                   WHERE gene_symbol LIKE %s {}
                  ORDER BY gene_symbol, organism_id, ensembl_release DESC
                """.format(org_addendum)
                cursor.execute(qry, ('%' + gene_symbol + '%',))

            for (id, ensembl_id, ensembl_version, ensembl_release, genbank_acc, organism_id,
                 molecule, start, stop, gene_symbol, product, biotype) in cursor:

                lc_gene_symbol = gene_symbol.lower()

                if organism_id in found[lc_gene_symbol]:
                    if found[lc_gene_symbol][organism_id]['ensembl_release'] == ensembl_release:
                        found[lc_gene_symbol][organism_id]['count'] += 1

                    continue
                else:
                    found[lc_gene_symbol][organism_id] = {'ensembl_release': ensembl_release,
                                                          'count': 1}

                gene = Gene(id=id, ensembl_id=ensembl_id, ensembl_version=ensembl_version,
                            ensembl_release=ensembl_release, genbank_acc=genbank_acc,
                            organism_id=organism_id, molecule=molecule, start=start,
                            stop=stop, gene_symbol=gene_symbol, product=product,
                            biotype=biotype)

                self.add_gene(gene)

        cursor.close()
        conn.close()


class GeneCart:
    def __init__(self, id=None, user_id=None, gctype=None, label=None, ldesc=None, organism_id=None,
                 genes=None, share_id=None, is_public=None, is_domain=None, date_added=None,
                 folder_id=None, folder_parent_id=None, folder_label=None):
        self.id = id
        self.user_id = user_id
        self.gctype = gctype
        self.label = label
        self.organism_id = organism_id
        self.ldesc = ldesc
        self.share_id = share_id
        self.is_public = is_public
        self.is_domain = is_domain
        self.date_added = date_added

        if not share_id:
            self.share_id = str(uuid.uuid4()).split('-')[0]

        # TODO: This should be a reference to a GeneCollection
        self.genes = list()
        self.num_genes = len(self.genes)

        # The are derived, populated by GeneCartCollection methods
        self.folder_id = folder_id
        self.folder_parent_id = folder_parent_id
        self.folder_label = folder_label

    def __repr__(self):
        return json.dumps(self.__dict__)

    def add_gene(self, gene):
        self.genes.append(gene)

    def get_genes(self):
        conn = Connection()
        cursor = conn.get_cursor()

        qry = "SELECT gene_symbol FROM gene_cart_member WHERE gene_cart_id = %s"
        cursor.execute(qry, (self.id,))

        self.genes = list()

        for row in cursor:
            self.genes.append(row[0])

        self.num_genes = len(self.genes)

        cursor.close()
        conn.close()

    def get_gene_counts(self):
        conn = Connection()
        cursor = conn.get_cursor()

        qry = "SELECT COUNT(*) FROM gene_cart_member WHERE gene_cart_id = %s"
        cursor.execute(qry, (self.id,))
        self.num_genes = cursor.fetchone()[0]

        cursor.close()
        conn.close()

    def remove(self):
        """
        Removes a gene cart and its dependencies from the database.  Requires object's
        'id' attribute to be defined.
        """
        if not self.id:
            raise Exception("Failed to delete a gene cart without an ID")

        # gene_cart_member entries are deleted by foreign key cascade
        conn = Connection()
        cursor = conn.get_cursor()

        sql = "DELETE FROM gene_cart WHERE id = %s"
        cursor.execute(sql, (self.id,))

        conn.commit()
        cursor.close()
        conn.close()

    def save(self):
        """
        Will perform a save or an update depending on whether the ID attribute is
        defined.
        """
        conn = Connection()
        cursor = conn.get_cursor()

        gcm_insert_qry = "INSERT INTO gene_cart_member (gene_cart_id, gene_symbol) VALUES (%s, %s)"

        if self.id is None:
            # ID is empty, this is a new one
            #  Insert the cart and then add the members
            gc_insert_qry = """
                            INSERT INTO gene_cart (user_id, label, organism_id, share_id, is_public, is_domain, gctype, ldesc)
                            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
            """

            cursor.execute(gc_insert_qry, (self.user_id, self.label, self.organism_id, self.share_id, self.is_public, self.is_domain, self.gctype, self.ldesc))
            self.id = cursor.lastrowid

            for gene in self.genes:
                cursor.execute(gcm_insert_qry, (self.id, gene.gene_symbol))

        else:
            # ID already populated
            #  Update cart properties, delete existing members, add current ones
            raise Exception("Called feature not yet implemented")

        cursor.close()
        conn.commit()
        conn.close()

    def save_change(self, attribute=None, value=None):
        """
        Update a cart attribute, both in the object and the relational database
        """
        if self.id is None:
            raise Exception("Error: no gene cart id. Cannot save change.")
        if attribute is None:
            raise Exception("Error: no attribute given. Cannot save change.")

        ## quick sanitization of attribute
        attribute = re.sub('[^a-zA-Z0-9_]', '_', attribute)
        setattr(self, attribute, value)

        conn = Connection()
        cursor = conn.get_cursor()

        save_sql = """
            UPDATE gene_cart
            SET {0} = %s
            WHERE id = %s
        """.format(attribute)
        cursor.execute(save_sql, (str(value), self.id))

        conn.commit()
        cursor.close()
        conn.close()

    def update_from_json(self, json_obj):
        """
        Takes an existing GeneCart object and updates/overrides any attributes provided by the
        passed data structure.

        Assumes json_obj is a parsed version of a string like this:
        {"label": "Test",
         "session_id": "e385305c-4387-433e-8b62-4bcf7c30ac52",
         "genes":
           [{"id": 2804, "gene_symbol": "Otor"},
            {"id": 2957, "gene_symbol": "Tgfbi"},
            {"id": 15772, "gene_symbol": "H19"},
            {"id": 168014, "gene_symbol": "Kctd12"}]
        }
        """
        if 'label' in json_obj:
            self.label = json_obj['label']

        if 'organism_id' in json_obj:
            self.organism_id = json_obj['organism_id']

        if 'session_id' in json_obj and not self.user_id:
            user_logged_in = get_user_from_session_id(json_obj['session_id'])
            self.user_id = user_logged_in.id

        if 'gctype' in json_obj:
            self.gctype = json_obj['gctype']

        if 'is_public' in json_obj:
            self.is_public = json_obj['is_public']

        if 'is_domain' in json_obj:
            self.is_domain = json_obj['is_domain']

        if 'genes' in json_obj:
            self.genes = list()

            for gene_dict in json_obj['genes']:
                gene = Gene(gene_symbol=gene_dict['gene_symbol'])
                self.add_gene(gene)


@dataclass
class GeneCartCollection:
    carts: List[GeneCart] = field(default_factory=list)

    # should gene-populating methods called to include gene members (lots of overhead)
    include_genes: bool = True

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        # Called when json modules attempts to serialize
        return self.__dict__

    def _row_to_cart_object(self, row):
        """
        Utility function so we don't have to repeat the SQL->Python object conversion
        """
        cart = GeneCart(
                    id=row['id'],
                    user_id=row['user_id'],
                    organism_id=row['organism_id'],
                    gctype=row['gctype'],
                    label=row['label'],
                    ldesc=row['ldesc'],
                    share_id=row['share_id'],
                    is_public=row['is_public'],
                    is_domain=row['is_domain'],
                    date_added=row['date_added'].isoformat(),
                    folder_id=row['folder_id'],
                    folder_parent_id=row['parent_id'],
                    folder_label=row['folder_label']
                )

        return cart

    def get_by_cart_ids(self, ids=[]):
        if not ids:
            raise Exception("No cart IDs provided to get_by_cart_ids")

        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = """
              SELECT gc.id, gc.user_id, gc.organism_id, gc.gctype, gc.label, gc.ldesc, gc.share_id,
                     gc.is_public, gc.is_domain, gc.date_added, f.id as folder_id, f.parent_id,
                     f.label as folder_label
                FROM gene_cart gc
                     LEFT JOIN folder_member fm ON fm.item_id=gc.id
                     LEFT JOIN folder f ON f.id=fm.folder_id
               WHERE gc.id = %s
                 AND (fm.item_type = 'genecart' or fm.item_type is NULL)
            ORDER BY gc.label
        """

        for id in ids:
            cursor.execute(qry, (id,))

            for row in cursor:
                cart = self._row_to_cart_object(row)

                if self.include_genes == True:
                    cart.get_genes()
                else:
                    cart.get_gene_counts()

                self.carts.append(cart)

        cursor.close()
        conn.close()
        return self.carts

    def get_by_share_ids(self, share_ids=[]):
        if not share_ids:
            raise Exception("No share_ids provided to get_by_share_ids")

        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = """
              SELECT gc.id, gc.user_id, gc.organism_id, gc.gctype, gc.label, gc.ldesc, gc.share_id,
                     gc.is_public, gc.is_domain, gc.date_added, f.id as folder_id, f.parent_id,
                     f.label as folder_label
                FROM gene_cart gc
                     LEFT JOIN folder_member fm ON fm.item_id=gc.id
                     LEFT JOIN folder f ON f.id=fm.folder_id
               WHERE gc.share_id = %s
                 AND (fm.item_type = 'genecart' or fm.item_type is NULL)
            ORDER BY gc.label
        """

        for share_id in share_ids:
            cursor.execute(qry, (share_id,))

            for row in cursor:
                cart = self._row_to_cart_object(row)

                if self.include_genes == True:
                    cart.get_genes()
                else:
                    cart.get_gene_counts()

                self.carts.append(cart)

        cursor.close()
        conn.close()
        return self.carts

    def get_by_user(self, user=None):
        if not user:
            raise Exception("User not provided to get_by_user")

        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = """
              SELECT gc.id, gc.user_id, gc.organism_id, gc.gctype, gc.label, gc.ldesc, gc.share_id,
                     gc.is_public, gc.is_domain, gc.date_added, f.id as folder_id, f.parent_id,
                     f.label as folder_label
                FROM gene_cart gc
                     LEFT JOIN folder_member fm ON fm.item_id=gc.id
                     LEFT JOIN folder f ON f.id=fm.folder_id
               WHERE gc.user_id = %s
                 AND (fm.item_type = 'genecart' or fm.item_type is NULL)
            ORDER BY gc.label
        """

        cursor.execute(qry, (user.id,))

        for row in cursor:
            cart = self._row_to_cart_object(row)

            if self.include_genes == True:
                cart.get_genes()
            else:
                cart.get_gene_counts()

            self.carts.append(cart)

        cursor.close()
        conn.close()
        return self.carts

    def get_by_user_groups(self, user=None):
        """
        Put here as it will be needed in the future. User groups not yet supported.
        """

        if not user:
            raise Exception("User not provided to get_by_user_groups")

        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = """
              SELECT gc.id, gc.user_id, gc.organism_id, gc.gctype, gc.label, gc.ldesc,
                     gc.share_id, gc.is_public, gc.is_domain, gc.date_added, f.id as folder_id,
                     f.parent_id, f.label as folder_label
                FROM ggroup g
                     JOIN user_group_membership ugm ON ugm.group_id=g.id
                     JOIN guser u ON u.id=ugm.user_id
                     JOIN gene_cart_group_membership gcgm ON gcgm.group_id=g.id
                     JOIN gene_cart gc ON gcgm.gene_cart_id=gc.id
                     LEFT JOIN folder_member fm ON fm.item_id=gc.id
                     LEFT JOIN folder f ON f.id=fm.folder_id
               WHERE u.id = %s
                 AND (fm.item_type = 'genecart' or fm.item_type is NULL)
            ORDER BY gc.label
        """
        cursor.execute(qry, (user.id,))

        for row in cursor:
            cart = self._row_to_cart_object(row)

            if self.include_genes == True:
                cart.get_genes()
            else:
                cart.get_gene_counts()

            self.carts.append(cart)

        cursor.close()
        conn.close()
        return self.carts

    def get_by_user_recent(self, user=None, n=None):
        if not user:
            raise Exception("No user provided to get_by_user_recent")

        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = """
              SELECT gc.id, gc.user_id, gc.organism_id, gc.gctype, gc.label, gc.ldesc, gc.share_id,
                     gc.is_public, gc.is_domain, gc.date_added, f.id as folder_id, f.parent_id,
                     f.label as folder_label
                FROM gene_cart gc
                     LEFT JOIN folder_member fm ON fm.item_id=gc.id
                     LEFT JOIN folder f ON f.id=fm.folder_id
               WHERE gc.user_id = %s
                 AND (fm.item_type = 'genecart' or fm.item_type is NULL)
            ORDER BY gc.date_added DESC
        """
        cursor.execute(qry, (user.id,))

        rows_returned = 0

        for row in cursor:
            cart = self._row_to_cart_object(row)

            if self.include_genes == True:
                cart.get_genes()
            else:
                cart.get_gene_counts()

            self.carts.append(cart)
            rows_returned += 1

            if n is not None and rows_returned >= n:
                break

        cursor.close()
        conn.close()
        return self.carts

    def get_domain(self):
        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = """
              SELECT gc.id, gc.user_id, gc.organism_id, gc.gctype, gc.label, gc.ldesc, gc.share_id,
                     gc.is_public, gc.is_domain, gc.date_added, f.id as folder_id, f.parent_id,
                     f.label as folder_label
                FROM gene_cart gc
                     LEFT JOIN folder_member fm ON fm.item_id=gc.id
                     LEFT JOIN folder f ON f.id=fm.folder_id
               WHERE gc.is_domain = 1
                 AND (fm.item_type = 'genecart' or fm.item_type is NULL)
            ORDER BY gc.label
        """
        cursor.execute(qry)

        for row in cursor:
            cart = self._row_to_cart_object(row)

            if self.include_genes == True:
                cart.get_genes()
            else:
                cart.get_gene_counts()

            self.carts.append(cart)

        cursor.close()
        conn.close()
        return self.carts

    def get_public(self):
        conn = Connection()
        cursor = conn.get_cursor(use_dict=True)

        qry = """
              SELECT gc.id, gc.user_id, gc.organism_id, gc.gctype, gc.label, gc.ldesc, gc.share_id,
                     gc.is_public, gc.is_domain, gc.date_added, f.id as folder_id, f.parent_id,
                     f.label as folder_label
                FROM gene_cart gc
                     LEFT JOIN folder_member fm ON fm.item_id=gc.id
                     LEFT JOIN folder f ON f.id=fm.folder_id
               WHERE gc.is_public = 1
                 AND (fm.item_type = 'genecart' or fm.item_type is NULL)
            ORDER BY gc.label
        """
        cursor.execute(qry)

        for row in cursor:
            cart = self._row_to_cart_object(row)

            if self.include_genes == True:
                cart.get_genes()
            else:
                cart.get_gene_counts()

            self.carts.append(cart)

        cursor.close()
        conn.close()
        return self.carts


class LayoutMember:
    def __init__(self, id=None, dataset_id=None, grid_position=None, mg_grid_position=None, start_col=None, mg_start_col=None, grid_width=None, mg_grid_width=None, start_row=None, mg_start_row=None, grid_height=None, mg_grid_height=None):
        self.id = id
        self.dataset_id = dataset_id
        self.grid_position = grid_position
        self.mg_grid_position = mg_grid_position
        self.start_col = start_col  # if not provided, fill into layout where it fits
        self.mg_start_col = mg_start_col
        self.grid_width = grid_width
        self.mg_grid_width = mg_grid_width
        self.start_row = start_row  # if not provided, fill into layout where it fits
        self.mg_start_row = mg_start_row
        self.grid_height = grid_height
        self.mg_grid_height = mg_grid_height

    def __repr__(self):
        return json.dumps(self.__dict__)

    def remove(self):
        """
        Deletes the current layout member from the database
        """
        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              DELETE FROM layout_members
              WHERE id = %s
        """
        cursor.execute(qry, (self.id,))

        cursor.close()
        conn.commit()
        conn.close()

    def save(self, layout=None):
        """
        Will perform a save or an update depending on whether the ID attribute is
        defined.  If ID is already present, it's assumed to be an existing row which
        needs to be updated.
        """
        if layout is None:
            raise Exception("You must pass a layout member object to save a layout member")

        conn = Connection()
        cursor = conn.get_cursor()

        if self.id is None:
            lm_insert_qry = """
            INSERT INTO layout_members (layout_id, dataset_id, grid_position, mg_grid_position, start_col, mg_start_col, grid_width, mg_grid_width, start_row, mg_start_row, grid_height, mg_grid_height)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
            """
            cursor.execute(lm_insert_qry, (layout.id, self.dataset_id, self.grid_position, self.mg_grid_position,
                                            self.start_col, self.mg_start_col, self.grid_width,
                                            self.mg_grid_width, self.start_row, self.mg_start_row,
                                            self.grid_height, self.mg_grid_height))
            self.id = cursor.lastrowid
        else:
            # ID already populated
            # Update Layout properties, delete existing members, add current ones
            lm_update_qry = """
            UPDATE layout_members
            SET dataset_id = %s, grid_position = %s, mg_grid_position = %s, start_col = %s, mg_start_col = %s, grid_width = %s, mg_grid_width = %s, start_row = %s, mg_start_row = %s, grid_height = %s, mg_grid_height = %s
            WHERE id = %s and layout_id = %s
            """
            cursor.execute(lm_update_qry, (self.dataset_id, self.grid_position, self.mg_grid_position,
                                           self.start_col, self.mg_start_col, self.grid_width,
                                           self.mg_grid_width, self.start_row, self.mg_start_row,
                                           self.grid_height, self.mg_grid_height, self.id, layout.id))

        cursor.close()
        conn.commit()
        conn.close()

class LayoutDisplay:
    def __init__(self, id=None, display_id=None, grid_position=None, start_col=None, grid_width=None, start_row=None, grid_height=None):
        self.id = id
        self.display_id = display_id
        self.dataset_id = None
        self.is_multigene = None  # True if multi, False if single
        self.grid_position = grid_position
        self.start_col = start_col  # if not provided, fill into layout where it fits
        self.grid_width = grid_width
        self.start_row = start_row  # if not provided, fill into layout where it fits
        self.grid_height = grid_height

    def __repr__(self):
        return json.dumps(self.__dict__)

    def remove(self):
        """
        Deletes the current layout member from the database
        """
        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              DELETE FROM layout_displays
              WHERE id = %s
        """
        cursor.execute(qry, (self.id,))

        cursor.close()
        conn.commit()
        conn.close()

    def get_dataset_id(self):
        conn = Connection()
        cursor = conn.get_cursor()

        qry = "SELECT dataset_id FROM dataset_display WHERE id = %s"
        cursor.execute(qry, (self.display_id,))

        for (dataset_id,) in cursor:
            self.dataset_id = dataset_id

        cursor.close()
        conn.close()

    def get_is_multigene(self):
        single_plot_types = ["bar", "line", "scatter", "tsne/umap_dynamic", "tsne_dynamic", "violin", "pca_static", "tsne_static", "tsne", "umap_static", "svg", "epiviz"]
        multi_plot_types = ["dotplot", "heatmap", "mg_violin", "quadrant", "volcano"]

        qry = """
            SELECT plot_type from dataset_display WHERE id = %s
            """
        conn = Connection()
        cursor = conn.get_cursor()
        cursor.execute(qry, (self.display_id,))
        (plot_type,) = cursor.fetchone()

        is_single = plot_type in single_plot_types
        is_multi = plot_type in multi_plot_types

        self.is_multigene = is_multi

        cursor.close()
        conn.close()


    def save(self, layout=None):
        """
        Will perform a save or an update depending on whether the ID attribute is
        defined.  If ID is already present, it's assumed to be an existing row which
        needs to be updated.
        """
        if layout is None:
            raise Exception("You must pass a layout member object to save a layout member")

        conn = Connection()
        cursor = conn.get_cursor()

        if self.id is None:
            lm_insert_qry = """
            INSERT INTO layout_displays (layout_id, display_id, grid_position, start_col, grid_width, start_row, grid_height)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
            """
            cursor.execute(lm_insert_qry, (layout.id, self.display_id, self.grid_position, self.start_col, self.grid_width, self.start_row, self.grid_height))
            self.id = cursor.lastrowid
        else:
            # ID already populated
            # Update Layout properties, delete existing members, add current ones
            lm_update_qry = """
            UPDATE layout_displays
            SET display_id = %s, grid_position = %s, start_col = %s, grid_width = %s, start_row = %s, grid_height = %s
            WHERE id = %s and layout_id = %s
            """
            cursor.execute(lm_update_qry, (self.display_id, self.grid_position, self.start_col, self.grid_width, self.start_row, self.grid_height, self.id, layout.id))

        cursor.close()
        conn.commit()
        conn.close()

class User:
    """
    Important note.  Because 'pass' is a reserved word in Python this field differs from the database
    table column name.
    """
    def __init__(self, id=None, user_name=None, email=None, institution=None, password=None,
                 updates_wanted=None, is_admin=None, default_org_id=None, is_curator=None,
                 help_id=None, layout_share_id=None, colorblind_mode=None):
        self.id = id
        self.user_name = user_name
        self.email = email
        self.institution = institution
        self.password = password
        self.colorblind_mode = colorblind_mode
        self.updates_wanted = updates_wanted
        self.is_admin = is_admin
        self.default_org_id = default_org_id
        self.is_curator = is_curator
        self.help_id = help_id
        self.layout_share_id = layout_share_id

        # This is a DatasetCollection and is NOT guaranteed to be all the user's datasets
        #  It depends on the options passed when calling self.datasets()
        self._datasets = None

        # This is a list of Layout objects
        self._layouts = None

    def __repr__(self):
        return json.dumps(self.__dict__)

    def datasets(self, has_h5ad=None, types=None):
        # populate it if we haven't already
        if self._datasets is None:
            self._datasets = DatasetCollection()
            self._datasets.get_owned_by_user(has_h5ad=has_h5ad, user=self, types=types)

        return self._datasets

