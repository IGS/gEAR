import datetime
import json
import os
import re
import sys
import uuid

from collections import defaultdict
from dataclasses import dataclass, field
from typing import List
from json import JSONEncoder

gear_lib_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(gear_lib_path)
import gear.db

# https://stackoverflow.com/a/35904211/1368079
this = sys.modules[__name__]

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

def _read_domain_short_label():
    json_conf = _read_site_domain_config()
    return json_conf['domain_short_display_label']

# For those functional differences we have depending on the site domain
#  Current values are gear, nemo.  Value kept in www/site_domain_prefs.json
this.domain_url = _read_domain_url()
this.domain_label = _read_domain_label()
this.domain_short_label = _read_domain_short_label()

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
         SELECT id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, ldesc, date_added,
                dtype, schematic_image, share_id, math_default, marked_for_removal, load_status,
                has_h5ad
           FROM dataset
          WHERE id = %s
    """
    cursor.execute(qry, (id, ))
    dataset = None

    for (id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, ldesc, date_added,
         dtype, schematic_image, share_id, math_default, marked_for_removal, load_status,
         has_h5ad) in cursor:
        dataset = Dataset(id=id, owner_id=owner_id, title=title, organism_id=organism_id,
                          pubmed_id=pubmed_id, geo_id=geo_id, is_public=is_public, ldesc=ldesc,
                          date_added=date_added, dtype=dtype, schematic_image=schematic_image,
                          share_id=share_id, math_default=math_default,
                          marked_for_removal=marked_for_removal, load_status=load_status,
                          has_h5ad=has_h5ad)

    if include_shape == '1':
        dataset.get_shape()

    cursor.close()
    conn.close()
    return dataset


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

def get_layout_by_id(layout_id):
    """
    Given a passed layout_id returns a Layout object with all attributes
    populated.  Returns None if no layout is found with that ID.
    """
    conn = Connection()
    cursor = conn.get_cursor()
    layout = None

    qry = """
          SELECT id, user_id, group_id, label, is_current, share_id
          FROM layout
          WHERE id = %s
    """
    cursor.execute(qry, (layout_id,))

    for (id, user_id, group_id, label, is_current, share_id) in cursor:
        layout = Layout(id=id, user_id=user_id, group_id=group_id,
                        label=label, is_current=is_current, share_id=share_id)
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

def get_user_by_id(user_id):
    """
       Given a user_id string this returns a User object with
       all attributes populated.  Returns None if no user with
       that ID is found.
    """
    conn = Connection()
    cursor = conn.get_cursor()

    qry = """
          SELECT g.id, g.user_name, g.email, g.institution, g.pass, g.updates_wanted,
                 g.is_admin, g.is_gear_curator, g.help_id
            FROM guser g
           WHERE g.id = %s
    """
    cursor.execute(qry, (user_id, ) )

    user = None
    for (id, user_name, email, institution, password, updates_wanted, is_admin, is_gear_curator, help_id) in cursor:
        user = User(id=id, user_name=user_name, email=email, institution=institution,
                    password=password, updates_wanted=updates_wanted, is_admin=is_admin,
                    is_gear_curator=is_gear_curator, help_id=help_id)
        break

    cursor.close()
    conn.close()
    return user

def get_user_from_session_id(session_id):
    """
       Given a session_id string this returns a User object with
       all attributes populated.  Returns None if no user for
       that session is known.
    """

    conn = Connection()
    cursor = conn.get_cursor()

    qry = """
          SELECT g.id, g.user_name, g.email, g.institution, g.pass, g.updates_wanted,
                 g.is_admin, g.is_gear_curator, g.help_id
            FROM guser g
                 JOIN user_session us ON g.id=us.user_id
           WHERE us.session_id = %s
    """
    cursor.execute(qry, (session_id, ) )

    user = None
    for (id, user_name, email, institution, password, updates_wanted, is_admin, is_gear_curator, help_id) in cursor:
        user = User(id=id, user_name=user_name, email=email, institution=institution,
                    password=password, updates_wanted=updates_wanted, is_admin=is_admin,
                    is_gear_curator=is_gear_curator, help_id=help_id)
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
        cursor.close()
        conn.close()
        return dict(
            id=id, dataset_id=dataset_id, user_id=user_id, label=label,
            plot_type=plot_type, plotly_config=json.loads(plotly_config)
        )
    except:
        cursor.close()
        conn.close()
        return None

def get_default_display(user_id, dataset_id):
    """Return user's display preference for given dataset."""
    conn = Connection()
    cursor = conn.get_cursor()
    qry = """
        SELECT display_id FROM dataset_preference
        where user_id = %s and dataset_id = %s
    """
    cursor.execute(qry, (user_id, dataset_id))
    try:
        (default_display_id,) = cursor.fetchone()
        cursor.close()
        conn.close()
        return default_display_id
    except:
        cursor.close()
        conn.close()
        # User has no display preference for this dataset
        return None

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
        return open(pipeline_file).read()

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
                    raise Exception("ERROR: base_path() called on Analysis object with no user_id attribute set.")

                # ./$user_id/$dataset_id/$analysis_id/$dataset_id.h5ad
                return "{0}/../www/analyses/by_user/{1}/{2}/{3}".format(this_dir, self.user_id, self.dataset_id, self.id)

            elif self.type == 'user_unsaved':
                if self.session_id is None:
                    raise Exception("ERROR: base_path() called on Analysis object with no session_id attribute set.")

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
        if current_user.is_gear_curator:
            self.vetting = 'gear'
            return 'gear'
        # Are the current user and dataset owner the same?
        elif self.user_id == current_user_id:
            self.vetting = 'owner'
            return 'owner'
        else:
            self.vetting = 'community'
            return 'community'


    def discover_type(self, current_user_id=None):
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
        ana = Analysis(
            id=jsn['id'], dataset_id=jsn['dataset_id'], label=jsn['label'], session_id=jsn['user_session_id'],
            user_id=None, type=jsn['type']
        )

        ## get the rest
        for k in jsn:
            if not hasattr(ana, k):
                # some were manually named, skip them
                if k not in ['user_session_id']:
                    setattr(ana, k, jsn[k])

        return ana

    def get_adata(self, backed=False, force_sparse=False):
        """
        Returns the anndata object for the current analysis.
        """
        # This goes against PEP8, but putting the import here should speed up the
        #  common case of most of the rest of this module where scanpy isn't needed
        import scanpy as sc

        # TODO: This could be ugly if we have to keep branching on read_h5ad options.  Clean it up.
        if backed:
            if force_sparse:
                return sc.read_h5ad(self.dataset_path(), backed='r', as_sparse="raw.X")
            else:
                return sc.read_h5ad(self.dataset_path(), backed='r')
        else:
            if force_sparse:
                return sc.read_h5ad(self.dataset_path(), as_sparse="raw.X")
            else:
                return sc.read_h5ad(self.dataset_path())


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
                            json_obj = json.loads(open(json_path).read())
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
        self.user_saved = self._scan_analysis_directory(ana, 'user_saved')
        self.user_unsaved = self._scan_analysis_directory(ana, 'user_unsaved')

class Connection:
    def __init__(self):
        self.mysql_cnx = gear.db.MySQLDB().connect()

    def commit(self):
        self.mysql_cnx.commit()

    def close(self):
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

class OrganismCollection:
    def __init__(self, organisms=None):
        self.organisms = [] if organisms is None else organisms

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

        return self.organisms
    
class Layout:
    def __init__(self, id=None, user_id=None, group_id=None, label=None,
                 is_current=None, share_id=None, members=None):
        self.id = id
        self.user_id = user_id
        self.group_id = group_id
        self.label = label
        self.is_current = is_current
        self.share_id = share_id
 
        # This should be a list of LayoutMember objects
        if not members:
            self.members = list()

        # handle defaults
        # TODO: If is_current = 1 we really need to reset all the other layouts by this user
        if not is_current:
            self.is_current = 0

        if not share_id:
            self.share_id = str(uuid.uuid4()).split('-')[0]

    def __repr__(self):
        return json.dumps(self.__dict__)

    def add_member(self, member):
        """
        Adds a LayoutMember to the database as part of this Layout
        """
        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              INSERT INTO layout_members (layout_id, dataset_id, grid_position, grid_width)
              VALUES (%s, %s, %s, %s)
        """
        cursor.execute(qry, (self.id, member.dataset_id, member.grid_position, member.grid_width))
        member.id = cursor.lastrowid
        self.members.append(member)

        cursor.close()
        conn.commit()

    def get_members(self):
        """
        Gets all members from the database and populates the 'members' attribute as a
        list of LayoutMember objects
        """
        conn = Connection()
        cursor = conn.get_cursor()

        self.members = list()

        qry = """
              SELECT lm.id, lm.dataset_id, lm.grid_position, lm.grid_width
                FROM layout_members lm
                     JOIN dataset d ON lm.dataset_id=d.id
               WHERE lm.layout_id = %s
                     AND d.marked_for_removal = 0
            ORDER BY lm.grid_position
        """
        cursor.execute(qry, (self.id,))

        for row in cursor:
            lm = LayoutMember(id=row[0], dataset_id=row[1], grid_position=row[2], grid_width=row[3])
            self.members.append(lm)

        cursor.close()

    def load(self):
        """
        If you only have the ID of layout and create an object from it, this method loads
        all the rest of the attributes, including layout members.
        """
        self.members = list()
        
        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              SELECT user_id, group_id, label, is_current, share_id
                FROM layout
               WHERE id = %s
        """
        cursor.execute(qry, (self.id,))

        for row in cursor:
            (self.user_id, self.group_id, self.label, self.is_current, self.share_id) = row

        self.get_members()
            
        cursor.close()
        conn.commit()

    def remove(self):
        """
        Deletes the current layout from the database, along with its layout members.
        """
        self.remove_all_members()

        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              DELETE FROM layout
              WHERE id = %s
        """
        cursor.execute(qry, (self.id,))

        cursor.close()
        conn.commit()

    def remove_all_members(self):
        """
        Deletes all members for this layout from the database in one query
        """
        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              DELETE FROM layout_members
              WHERE layout_id = %s
        """
        cursor.execute(qry, (self.id,))
        
        cursor.close()
        conn.commit()

        self.members = []

    def remove_member_by_dataset_id(self, dataset_id):
        """
        Rather than the layout_member.id, this is a utility function to delete
        a member based on its dataset ID.  In the future it will be possible for a
        profile to have two representations of the same dataset in the profile, and
        this will have to go.  For now we do a check to protect for this just in case.
        """
        found_member_dataset_ids = list()
        for lm in self.members:
            if lm.dataset_id in found_member_dataset_ids:
                raise Exception("ERROR: Found two datasets with same ID as part of the same layout.  Not safe to remove based on dataset_id alone.")
            else:
                found_member_dataset_ids.append(lm.dataset_id)

        conn = Connection()
        cursor = conn.get_cursor()

        qry = """
              DELETE FROM layout_members
              WHERE dataset_id = %s
                AND layout_id = %s
        """
        cursor.execute(qry, (dataset_id, self.id))

        # make sure this member is removed from our internal list too
        self.members = [i for i in self.members if i.dataset_id != dataset_id]

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
            INSERT INTO layout (user_id, group_id, label, is_current, share_id)
            VALUES (%s, %s, %s, %s, %s)
            """
            cursor.execute(layout_insert_qry, (self.user_id, self.group_id, self.label,
                                               self.is_current, self.share_id))
            self.id = cursor.lastrowid
        else:
            # ID already populated
            # Update Layout properties, delete existing members, add current ones
            raise Exception("Layout.save() not yet implemented for update mode")

        cursor.close()
        conn.commit()

    # TODO: Need a function to take a DatasetCollection and populate
    #  information on it within a layout

        
@dataclass
class Dataset:
    id: str
    owner_id: int = None
    title: str = None
    organism_id: int = None
    pubmed_id: str = None
    geo_id: str = None
    is_public: int = None
    ldesc: str = None
    date_added: datetime.datetime = None
    dtype: str = None
    schematic_image: str = None
    share_id: str = None
    math_default: str = None
    marked_for_removal: int = None
    load_status: str = None
    has_h5ad: int = None
    platform_id: str = None
    instrument_model: str = None
    library_selection: str = None
    library_source: str = None
    library_strategy: str = None
    contact_email: str = None
    contact_institute: str = None
    contact_name: str = None
    annotation_source: str = None
    plot_default: str = None
    annotation_release: int = None
    # derived, here for convenience
    gene_count: int = None
    obs_count: int = None
    tags: List[str] = field(default_factory=list)
    layouts: List[Layout] = field(default_factory=list)

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _serialize_json(self):
        return self.__dict__

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
        if self.has_h5ad:
            if session_id is None:
                h5ad_file_path = "{0}/../www/datasets/{1}.h5ad".format(
                    os.path.dirname(os.path.abspath(__file__)), self.id)
            else:
                h5ad_file_path = "{0}/{1}/{2}.h5ad".format(this.analysis_base_dir, session_id, self.id)

            return h5ad_file_path
        else:
            ## all other types are in the same place
            tab_file_path = "{0}/../www/datasets_uploaded/{1}.tab".format(
                os.path.dirname(os.path.abspath(__file__)), self.id)

            return tab_file_path

    def get_layouts(self, user=None):
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
                      SELECT l.id, l.user_id, l.group_id, l.label, l.is_current, l.share_id
                        FROM layout l
                             JOIN layout_members lm ON lm.layout_id=l.id
                       WHERE lm.dataset_id = %s
                             AND (user_id = 0 OR user_id = %s)
                    ORDER BY l.label
                """
                cursor.execute(qry, (self.id, user.id))
            else:
                qry = """
                      SELECT l.id, l.user_id, l.group_id, l.label, l.is_current, l.share_id
                        FROM layout l
                             JOIN layout_members lm ON lm.layout_id=l.id
                       WHERE lm.dataset_id = %s
                             AND user_id = 0
                    ORDER BY l.label
                """
                cursor.execute(qry, (self.id,))

            for row in cursor:
                l = Layout(id=row[0], user_id=row[1], group_id=row[2], label=row[3],
                           is_current=row[4], share_id=row[5])
                self.layouts.append(l)
                
            cursor.close()
                
        return self.layouts
        
        
    def get_shape(self, session_id=None):
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

            import scanpy as sc
            sc.settings.verbosity = 0
            adata = sc.read_h5ad(h5ad_file_path)
            (n_obs, n_vars) = adata.shape

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

    def get_by_dataset_ids(self, ids=None):
        conn = Connection()
        cursor = conn.get_cursor()

        # I tried but couldn't make this work in a single query with IN () instead.
        qry = """
           SELECT d.id, d.title, o.label, d.pubmed_id, d.geo_id, d.is_public, d.ldesc,
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
        GROUP BY d.id, d.title, o.label, d.pubmed_id, d.geo_id, d.is_public, d.ldesc,
                   d.dtype, u.id, u.user_name, d.schematic_image, d.share_id,
                   d.math_default, d.marked_for_removal, d.date_added, d.load_status,
                   d.annotation_source, d.annotation_release, d.organism_id
        """
        for id in ids:
            cursor.execute(qry, (id,))

            for row in cursor:
                # skip datasets marked_for_removal
                if row[13] == 1:
                    continue
                else:
                    if row[5] == 1:
                        access_level = 'Public'
                    else:
                        access_level = 'Private'

                    date_added = row[14].isoformat()

                    if row[16] == 'NULL':
                        tag_list = None
                    else:
                        tag_list = row[16].replace(',', ', ')

                    dataset = Dataset(id=row[0],
                                      title=row[1],
                                      organism_id=row[19],
                                      pubmed_id=row[3],
                                      geo_id=row[4],
                                      is_public=row[5],
                                      ldesc=row[6],
                                      dtype=row[7],
                                      owner_id=row[8],
                                      schematic_image=row[10],
                                      share_id=row[11],
                                      math_default=row[12],
                                      date_added=date_added,
                                      load_status=row[15],
                                      annotation_source=row[17],
                                      annotation_release=row[18]
                    )
                        
                    # Add supplemental attributes this method created previously
                    dataset.organism = row[2]
                    dataset.tags = row[16].split(',')
                    dataset.access = 'access_level'
                    dataset.user_name = row[9]
                    
                    #  TODO: These all need to be tracked through the code and removed
                    dataset.dataset_id = dataset.id
                    dataset.user_id = dataset.owner_id
                    dataset.math_format = dataset.math_default

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
              SELECT id, owner_id, title, organism_id, pubmed_id, geo_id, is_public, ldesc,
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
                              geo_id=r['geo_id'], is_public=r['is_public'], ldesc=r['ldesc'],
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
              SELECT d.id, d.owner_id, d.title, organism_id, pubmed_id, geo_id, is_public, ldesc,
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
                              geo_id=r['geo_id'], is_public=r['is_public'], ldesc=r['ldesc'],
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
              SELECT d.id, d.owner_id, d.title, organism_id, pubmed_id, geo_id, is_public, ldesc,
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
                              geo_id=r['geo_id'], is_public=r['is_public'], ldesc=r['ldesc'],
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
    id: int = None
    ensembl_id: str = None
    ensembl_version: str = None
    ensembl_release: int = None
    genbank_acc: str = None
    organism_id: int = None
    molecule: str = None
    start: int = None
    stop: int = None
    gene_symbol: str = None
    product: str = None
    biotype: str = None

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

        if this.domain_label == 'gear':
            # Add one for the SHIELD
            self.dbxrefs.append({'source': 'SHIELD', 'identifier': self.gene_symbol, 'url': None})

            # Hack to add DVD links until they have SSL enabled and we can use their API instead
            dvd_genes = ["ACTG1","ADCY1","ADGRV1","GPR98","AIFM1","ALMS1","ATP2B2","ATP6V1B1","BDP1","BSND","C10orf2","CABP2","CACNA1D","CCDC50","CD164","CDC14A","CDH23","CEACAM16","CIB2","CISD2","CLDN14","CLIC5","CLPP","CLRN1","COCH","COL11A1","COL11A2","COL2A1","COL4A3","COL4A4","COL4A5","COL4A6","COL9A1","COL9A2","CRYM","DCDC2","DFNA5","DFNB31 (WHRN)","DFNB59 (PJVK)","DIABLO","DIAPH1","DIAPH3","DSPP","EDN3","EDNRB","ELMOD3","EPS8","EPS8L2","ESPN","ESRRB","EYA1","EYA4","FAM65B","FGF3","FGFR1","FGFR2","FOXI1","GATA3","GIPC3","GJB2","GJB3","GJB6","GPSM2","GRHL2","GRXCR1","GRXCR2","HARS2","HGF","HOMER2","HSD17B4","ILDR1","KARS","KCNE1","KCNJ10","KCNQ1","KCNQ4","KITLG","LARS2","LHFPL5","LOXHD1","LOXL3","LRTOMT","MARVELD2","MCM2","MET","MIR96","MITF","MSRB3","MT-RNR1","MT-TL1","MT-TS1","MYH14","MYH9","MYO15A","MYO3A","MYO6","MYO7A","NARS2","NLRP3","OPA1","OSBPL2","OTOA","OTOF","OTOG","OTOGL","P2RX2","PAX3","PCDH15","PDZD7","PEX1","PEX6","PNPT1","POLR1C","POLR1D","POU3F4","POU4F3","PRPS1","PTPRQ","RDX","ROR1","S1PR2","SERPINB6","SIX1","SIX5","SLC17A8","SLC22A4","SLC26A4","SLC26A5","SLITRK6","SMPX","SNAI2","SOX10","STRC","SYNE4","TBC1D24","TBX1","TCOF1","TECTA","TECTB","TIMM8A","TJP2","TMC1","TMEM132E","TMIE","TMPRSS3","TNC","TPRN","TRIOBP","TSPEAR","USH1C","USH1G","USH2A","WFS1"]
            if self.gene_symbol.upper() in dvd_genes:
                dvd_url = 'http://deafnessvariationdatabase.org/gene/' + self.gene_symbol.upper()
                self.dbxrefs.append({'source': 'DVD', 'identifier': self.gene_symbol, 'url': dvd_url, 'title': 'Deafness Variation Database'})

        elif this.domain_label == 'nemo':
            # Add one for UCSC Cell Browser, BrainSpan
            self.dbxrefs.append({'source': 'UCSC Cell Browser', 'identifier': self.gene_symbol.upper(), 'url': None})
            self.dbxrefs.append({'source': 'BrainSpan', 'identifier': self.gene_symbol.upper(), 'url': None})
            self.dbxrefs.append({'source': 'Cortecon', 'identifier': self.gene_symbol.upper(), 'url': None})

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
            elif source == 'UCSC Cell Browser':
                dbxref['url'] = "https://cells.ucsc.edu/?ds=cortex-dev&gene={0}".format(dbxref['identifier'])
            elif source == 'BrainSpan':
                dbxref['url'] = "http://www.brainspan.org/rnaseq/searches?exact_match=true&search_term={0}&search_type=gene".format(dbxref['identifier'])
            elif source == 'Cortecon':
                dbxref['url'] = "http://cortecon.neuralsci.org/index.php?genestring={0}&cort_mode=genesearch".format(dbxref['identifier'])
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

    def get_by_gene_symbol(self, gene_symbol=None, exact=None):
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

        for gene_symbol in gene_symbols:
            if len(gene_symbol) == 0:
                continue

            # this comes from javascript as a string true/false
            if exact in [True, 'true', 1, "1"]:
                qry = """
                  SELECT id, ensembl_id, ensembl_version, ensembl_release, genbank_acc, organism_id,
                         molecule, start, stop, gene_symbol, product, biotype
                    FROM gene
                   WHERE gene_symbol = %s
                  ORDER BY gene_symbol, organism_id, ensembl_release DESC
                """
                cursor.execute(qry, (gene_symbol,))
            else:
                qry = """
                  SELECT id, ensembl_id, ensembl_version, ensembl_release, genbank_acc, organism_id,
                         molecule, start, stop, gene_symbol, product, biotype
                    FROM gene
                   WHERE gene_symbol LIKE %s
                  ORDER BY gene_symbol, organism_id, ensembl_release DESC
                """
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
                 genes=None, share_id=None, is_public=None, date_added=None):
        self.id = id
        self.user_id = user_id
        self.gctype = gctype
        self.label = label
        self.organism_id = organism_id
        self.ldesc = ldesc
        self.share_id = share_id
        self.is_public = is_public
        self.date_added = date_added

        # TODO: This should be a reference to a GeneCollection
        if not genes:
            self.get_genes()

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
            gc_insert_qry = "INSERT INTO gene_cart (user_id, label, organism_id, share_id, is_public) VALUES (%s, %s, %s, %s, %s)"

            cursor.execute(gc_insert_qry, (self.user_id, self.label, self.organism_id, self.share_id, self.is_public))
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

        if 'genes' in json_obj:
            self.genes = list()

            for gene_dict in json_obj['genes']:
                gene = Gene(gene_symbol=gene_dict['gene_symbol'])
                self.add_gene(gene)


class LayoutMember:
    def __init__(self, id=None, dataset_id=None, grid_position=None, grid_width=None):
        self.id = id
        self.dataset_id = dataset_id
        self.grid_position = grid_position
        self.grid_width = grid_width

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
            INSERT INTO layout_members (layout_id, dataset_id, grid_position, grid_width)
            VALUES (%s, %s, %s, %s)
            """
            cursor.execute(lm_insert_qry, (layout.id, self.dataset_id, self.grid_position,
                                           self.grid_width))
            self.id = cursor.lastrowid
        else:
            # ID already populated
            # Update Layout properties, delete existing members, add current ones
            raise Exception("LayoutMember.save() not yet implemented for update mode")

        cursor.close()
        conn.commit()
        conn.close()


class User:
    """
    Important note.  Because 'pass' is a reserved word in Python this field differs from the database
    table column name.
    """
    def __init__(self, id=None, user_name=None, email=None, institution=None, password=None, updates_wanted=None,
                 is_admin=None, is_gear_curator=None, help_id=None):
        self.id = id
        self.user_name = user_name
        self.email = email
        self.institution = institution
        self.password = password
        self.updates_wanted = updates_wanted
        self.is_admin = is_admin
        self.is_gear_curator = is_gear_curator
        self.help_id = help_id

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

