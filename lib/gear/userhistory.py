import geardb
import re
import urllib

class UserHistory:
    def __init__(self, connection=False):
        """
        If you pass a geardb.Connection object it will be used to perform the query
        and committed but NOT closed. Otherwise a connection will be created automatically.
        """
        self.connection = connection

        if not self.connection:
            self.connection = geardb.Connection()

    def add_record(self, user_id=None, entry_category=None, label=None, **kwargs):
        """
        Adds a history record for a user action to the DB. Arguments needed for all types:

            user_id, entry_category, label

        Individual additional arguments needed depending on category:
        
            'dataset_search' -> search_terms (an array of terms)
            'gene_search' -> gene_symbol (can be a string with multiple), layout_share_id / dataset_share_id
            'layout_added' -> layout_share_id
            'multigene_search' -> gene_symbol (can be a string with multiple), layout_share_id / dataset_share_id
            'projection_run' -> patterns, algo, gene_cart, multi, layout_share_id

        """
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
                # looks like: https://nemoanalytics.org/p?p=p&multi=0&l=4e8f6c00&c=00be4b21&algo=pca
                if 'layout_share_id' in kwargs:
                    url = "/p?p=p&l={0}".format(kwargs['layout_share_id'])
                else:
                    raise Exception("ERROR: If recording a layout_added category, 'layout_share_id' must be passsed")
                
            case _ :
                raise Exception("ERROR: Invalid entry_category when calling UserHistory.add_record()")
        
        qry = """
              INSERT INTO user_history (user_id, entry_category, label, url)
              VALUES (%s, %s, %s, %s)
        """
        cursor = self.connection.get_cursor()
        cursor.execute(qry, (user_id, entry_category, label, url))
        self.connection.commit()

        
