import os, sys
import tarfile

class FileType(object):
    # filetypes = []
    def __init__(self, filetypes=None, adata=None):
        self.filetypes = filetypes
        self.adata = adata

        if self.filetypes is None:
            self.filetypes = list()

        #TODO: .adata added via ExcelUploader child-class
        if self.adata is None:
            pass

class DatasetUploader:

    def get_by_filetype(self, filetype=None, filepath=None):
        """
        This factory nests the dataset filetype classes. Preventing them from being directly called.

            dataset = DatasetUploader.get_by_filetype('xlsx')
            dataset.read_file(filepath)

        Creates an Excel class object which now can be used for process and upload the dataset.

        Follows example:
        http://python-3-patterns-idioms-test.readthedocs.io/en/latest/Factory.html#preventing-direct-creation
        """

        if filetype == "xlsx" or filetype == "xls":
            import gear.exceluploader as exceluploader
            return exceluploader.ExcelUploader()

        if filetype == "tar" or filetype == "gz":
            if filepath is None:
                raise Exception("filepath is a required argument for tarball uploading")
            else:
                tartype = self.tarball_type(filepath)
                if tartype == 'mex':
                    import gear.mexuploader as mexuploader
                    return mexuploader.MexUploader()
                elif tartype == 'threetab':
                    import gear.threetabuploader as threetabuploader
                    return threetabuploader.ThreeTabUploader()
                elif tartype == 'visium':
                    import gear.spatialuploader as spatialuploader
                    return spatialuploader.VisiumUploader()
                elif tartype == 'visiumhd':
                    import gear.spatialuploader as spatialuploader
                    return spatialuploader.VisiumHDUploader()
        if filetype == "h5ad":
            import gear.h5aduploader as h5aduploader
            return h5aduploader.H5adUploader()
        assert 0, "Do not recognize file type given: " + filetype

    def tarball_type(self, file=None):
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
        t = tarfile.open(file, 'r')
        filenames = t.getnames()
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

        # ? subject to change
        if "spatial" in filestr and "tissue_positions_list.csv" in filestr and "filtered_feature_bc_matrix.h5" in filestr:
            return 'visium'
        if "spatial" in filestr and "tissue_positions.parquet" in filestr and "filtered_feature_bc_matrix.h5" in filestr:
            return 'visiumhd'
        return None
