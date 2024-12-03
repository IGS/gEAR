import os, sys

import pandas as pd
import scanpy as sc
from scipy import sparse
import anndata
import json

# This has a huge dependency stack of libraries. Occasionally, one of them has methods
#  which prints debugging information on STDOUT, killing this CGI.  So here we redirect
#  STDOUT until we need it.
print('Content-Type: application/json\n\n', flush=True)
original_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')

class DatasetUploader:
    """
    This factory nests the dataset filetype classes. Preventing them from being directly called.

        dataset = DatasetUploader.get_by_filetype('xlsx')
        dataset.read_file(filepath)

    Creates an Excel class object which now can be used for process and upload the dataset.

    Follows example:
    http://python-3-patterns-idioms-test.readthedocs.io/en/latest/Factory.html#preventing-direct-creation
    """

    def __init__(self, share_uid=None, session_id=None, dataset_format=None, 
                 status_json_file=None, upload_dir=None, filetypes=None, adata=None):
        self.share_uid = share_uid
        self.session_id = session_id
        self.dataset_format = dataset_format
        self.status_json_file = status_json_file
        self.upload_dir = upload_dir

        self.status = {
            "success": None,
            "process_id": None,
            "status": "uploaded",
            "message": "",
            "progress": 0
        }

        self.filetypes = filetypes
        self.adata = adata

        if self.filetypes is None:
            self.filetypes = list()

        # TODO: .adata added via ExcelUploader child-class
        if self.adata is None:
            pass

    @classmethod
    def get_by_filetype(cls, **kwargs):
        filetype = kwargs.get('dataset_format', None)
        upload_dir = kwargs.get('upload_dir', None)
        share_uid = kwargs.get('share_uid', None)

        if filetype == "excel":
            import gear.exceluploader as exceluploader
            return exceluploader.ExcelUploader(**kwargs)

        if filetype == "mex_3tab":
            filepath = os.path.join(upload_dir, f"{share_uid}.tar.gz")
            uploaded_file_names = None
            file_found = False

            if os.path.exists(filepath):
                file_found = True
                compression_format = 'tarball'
                import tarfile

                t = tarfile.open(filepath, 'r')
                uploaded_file_names = t.getnames()
            else:
                filepath = os.path.join(upload_dir, f"{share_uid}.zip")

                if os.path.exists(filepath):
                    file_found = True
                    compression_format = 'zip'
                    import zipfile

                    z = zipfile.ZipFile(filepath, 'r')
                    uploaded_file_names = z.namelist()

            tartype = cls.tarball_type(uploaded_file_names)
            if tartype == 'mex':
                import gear.mexuploader as mexuploader
                uploader = mexuploader.MexUploader(**kwargs)
            elif tartype == 'threetab':
                import gear.threetabuploader as threetabuploader
                uploader = threetabuploader.ThreeTabUploader(**kwargs)

            if file_found == False:
                uploader.write_status(label='error', message="No tarball or zip file found.")
                uploader.dump_and_exit()
            else:
                return uploader
                
        if filetype == "rdata":
            import gear.rdatauploader as rdatauploader
            return rdatauploader.RdataUploader(**kwargs)
        
        if filetype == "h5ad":
            import gear.h5aduploader as h5aduploader
            return h5aduploader.H5adUploader(**kwargs)
        assert 0, "Do not recognize file type given: " + filetype

    @classmethod
    def tarball_type(cls, filenames):
        """
        Since we allow users to upload tarballs (.tar or .tar.gz) we don't know
        from the filename alone if this is MEX or ThreeTab.  This method reads the contents
        of the tarball and returns either 'mex' or 'threetab' based on what it finds. The basenames variable
        contains the filenames in case the files are within a forlder inside the tarball.

        mex:
        matrix.mtx
        barcodes.tsv
        genes.tsv

        threetab:
        expression.tab
        genes.tab
        observations.tab

        None is returned if neither of these is true
        
        Added NEMO file format functionality.
        DataMTX.tab -> expression.tab
        COLmeta.tab -> observations.tab
        ROWmeta.tab -> genes.tab
        """
        filestr=''.join(filenames)
        basenames=[]
        for f_path in filenames:
            fname=os.path.basename(f_path)
            basenames.append(fname)
            
        if 'expression.tab' in basenames and 'genes.tab' in basenames and 'observations.tab' in basenames:
            return 'threetab'
        
        if 'matrix.mtx' in basenames and 'barcodes.tsv' in basenames and 'genes.tsv' in basenames:
            return 'mex'
        
        if 'DataMTX.tab' in filestr and 'COLmeta.tab' in filestr and 'ROWmeta.tab' in filestr:
            return 'threetab'

        return None

    def write_status(self, success=None, label=None, message=None):
        if success is not None:
            if success:
                self.status['success'] = 1
            else:
                self.status['success'] = 0

        if label is not None:
            self.status['status'] = label

        if message is not None:
            self.status['message'] = message

        with open(self.status_json_file, 'w') as f:
            f.write(json.dumps(self.status))


    def dump_and_exit(self):
        sys.stdout = original_stdout

        response = { success: self.status['success'], message: self.status['message'] }
        print(json.dumps(response), flush=True)

        sys.exit(0)