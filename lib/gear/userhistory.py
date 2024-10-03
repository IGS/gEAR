import geardb
import re
import sys
import urllib
from json import JSONEncoder

from gear.serverconfig import ServerConfig

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

class UserHistory:
    """
    Represents a single entry in a user's activity history
    """
    def __init__(self, connection=False, user_id=None, entry_category=None,
                 entry_date=None, label=None, url=None):
        """
        If you pass a geardb.Connection object it will be used to perform the query
        and committed but NOT closed. Otherwise a connection will be created automatically.
        """
        self.connection = connection

        if not self.connection:
            self.connection = geardb.Connection()

        self.user_id = user_id
        self.entry_category = entry_category
        self.entry_date = entry_date
        self.label = label
        self.url = url

        server_config = ServerConfig().parse()
        self.user_history_enabled = server_config.getboolean('history', 'enable_user_history')
        self.history_count = server_config.getint('history', 'history_count')

    def _serialize_json(self):
        # Called when json modules attempts to serialize
        return self.__dict__

    def add_record(self, user_id=None, entry_category=None, label=None, **kwargs):
        """
        Returns False unless user_history is enabled in the server config.

        Adds a history record for a user action to the DB. Arguments needed for all types:

            user_id, entry_category, label

        Individual additional arguments needed depending on category:

            'dataset_search' -> search_terms (an array of terms)
            'gene_search' -> gene_symbol (can be a string with multiple), layout_share_id / dataset_share_id
            'layout_added' -> layout_share_id
            'multigene_search' -> gene_symbol (can be a string with multiple), layout_share_id / dataset_share_id
            'projection_run' -> patterns, algo, gene_cart, multi, layout_share_id

        """
        if not self.user_history_enabled:
            return False

        match entry_category:
            case 'dataset_search':
                if 'search_terms' in kwargs:
                    url = "/p?p=de&ss={}".format(urllib.parse.quote(str(" ".join(kwargs['search_terms']))))
                else:
                    raise Exception("ERROR: If recording a dataset_search category, 'search_terms' must be passsed")

            case 'gene_cart_added':
                if 'gene_cart_share_id' in kwargs:
                    url = "/p?p=gcm&s={0}".format(kwargs['gene_cart_share_id'])
                else:
                    raise Exception("ERROR: If recording a gene_cart_added category, 'gene_cart_share_id' must be passsed")

            case 'gene_cart_search':
                if 'search_terms' in kwargs:
                    url = "/p?p=gcm&ss={}".format(urllib.parse.quote(str(" ".join(kwargs['search_terms']))))
                else:
                    raise Exception("ERROR: If recording a dataset_search category, 'search_terms' must be passsed")

            case 'gene_search':
                if 'layout_share_id' in kwargs:
                    url = "/p?l={0}".format(kwargs['layout_share_id'])
                elif 'dataset_share_id' in kwargs:
                    url = "/p?s={0}".format(kwargs['dataset_share_id'])
                else:
                    raise Exception("ERROR: If recording a gene_search category, 'layout_share_id' or 'dataset_share_id' must be passsed")

                if 'gene_symbol' not in kwargs:
                    raise Exception("ERROR: If recording a gene_search category, 'gene_symbol' must be passed")

                gene_string = re.sub("[\, ]+", ",", kwargs['gene_symbol'])

                url += "&g={0}".format(gene_string)

            case 'layout_added':
                if 'layout_share_id' in kwargs:
                    url = "/p?l={0}".format(kwargs['layout_share_id'])
                else:
                    raise Exception("ERROR: If recording a layout_added category, 'layout_share_id' must be passsed")

            case 'multigene_search':
                if 'layout_share_id' in kwargs:
                    url = "/p?l={0}".format(kwargs['layout_share_id'])
                elif 'dataset_share_id' in kwargs:
                    url = "/p?s={0}".format(kwargs['dataset_share_id'])
                else:
                    raise Exception("ERROR: If recording a multigene_search category, 'layout_share_id' or 'dataset_share_id' must be passsed")

                if 'gene_symbol' not in kwargs:
                    raise Exception("ERROR: If recording a multigene_search category, 'gene_symbol' must be passed")

                gene_string = re.sub("[\, ]+", ",", kwargs['gene_symbol'])

                url += "&g={0}&multi=1&gsem=1".format(gene_string)

            case 'projection_run':
                # looks like: https://nemoanalytics.org/p?p=p&multi=0&l=4e8f6c00& c=00be4b21 &algo=pca
                if 'layout_share_id' in kwargs:
                    url = "/p?p=p&l={0}".format(kwargs['layout_share_id'])
                else:
                    raise Exception("ERROR: If recording a projection_run category, 'layout_share_id' must be passsed")

                if 'multi' in kwargs:
                    url += "&multi={0}".format(kwargs['multi'])
                else:
                    raise Exception("ERROR: If recording a projection_run category, 'multi' must be passsed")

                if 'gene_cart' in kwargs:
                    url += "&c={0}".format(kwargs['gene_cart'])
                else:
                    raise Exception("ERROR: If recording a projection_run category, 'gene_cart' must be passsed")

            case _ :
                raise Exception("ERROR: Invalid entry_category when calling UserHistory.add_record()")

        qry = """
              INSERT INTO user_history (user_id, entry_category, label, url)
              VALUES (%s, %s, %s, %s)
        """
        cursor = self.connection.get_cursor()
        cursor.execute(qry, (user_id, entry_category, label, url))
        self.connection.commit()

    def get_latest_entries(self, entry_category=None, entry_count=5, **kwargs):
        """
        Returns a list of UserHistory objects, usually for the purpose of populating something
        like a 'most recent activities' table.

        Can filter on category and/or number of entries to be returned.
        """
        cursor = self.connection.get_cursor()
        qry_args = [self.user_id]
        qry = """
                 SELECT entry_date, entry_category, label, url
                   FROM user_history
                  WHERE user_id = %s
                  ORDER BY entry_date DESC
        """

        if entry_category:
            qry += " AND entry_category = %s"
            qry_args.append(self.entry_category)

        qry += f" LIMIT {entry_count}"

        cursor.execute(qry, qry_args)
        entries = list()

        for row in cursor:
            entry = UserHistory(user_id=self.user_id, entry_category=row[1],
                                entry_date=row[0], label=row[2], url=row[3]
            )
            entries.append(entry)

        return entries












