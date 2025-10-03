import json
import sys
import typing
from pathlib import Path

# append parent directory to path
sys.path.append(str(Path(__file__).resolve().parents[1]))   # lib

from geardb import Dataset, get_user_by_id, get_user_id_from_session_id

if typing.TYPE_CHECKING:
    from anndata import AnnData
    from spatialdata import SpatialData

this_dir = Path(__file__).resolve().parent  # lib/gear
root_dir = this_dir.parents[1]

# You can use the analysis functions to create an Analysis|SpatialAnalysis object in a high-level way
# without needing to know the details of the class or its path conventions

# Creating the Analysis/SpatialAnalysis objects directly is useful if you are using theoretical paths instead
# of existing ones, such as copying an analysis. However the paths are pre-determined.

# If you just need to get the AnnData or SpatialData object for a custom path (i.e. the upload directory),
# you can use the adapters directly (H5adAdapter or ZarrAdapter)


def get_analysis(analysis_data: dict | None, dataset_id: str, session_id: str | None, is_spatial: bool = False) -> "SpatialAnalysis | Analysis":
    """
    Retrieves an analysis object (either SpatialAnalysis or Analysis) based on the provided analysis data, dataset ID, and session ID.

    If analysis data is provided, constructs the appropriate analysis object using the given information and verifies the existence of the associated h5ad file.
    If analysis data is not provided, returns the primary analysis for the specified dataset.

    Args:
        analysis_data (dict | None): Dictionary containing analysis metadata, or None to fetch the primary analysis.
        dataset_id (str): The unique identifier for the dataset.
        session_id (str): The session identifier for the current user.
        is_spatial (bool, optional): Whether to return a SpatialAnalysis object. Defaults to False.

    Returns:
        SpatialAnalysis | Analysis: The constructed or retrieved analysis object.

    Raises:
        FileNotFoundError: If the h5ad file for the specified analysis does not exist.
    """

    # If an analysis is posted we want to read from its h5ad
    if analysis_data:
        user_id = get_user_id_from_session_id(session_id)

        if is_spatial:
            ana = SpatialAnalysis(
                id=analysis_data["id"],
                dataset_id=dataset_id,
                session_id=session_id,
                user_id=user_id,
                type=analysis_data.get("type", None),
            )
        else:
            ana = Analysis(
                id=analysis_data["id"],
                dataset_id=dataset_id,
                session_id=session_id,
                user_id=user_id,
                type=analysis_data.get("type", None),
            )

        if ana.type is None:
            ana.discover_type()

        # Check that the h5ad file exists
        if not ana.dataset_path.exists():
            filetype = "zarr" if is_spatial else "h5ad"
            raise FileNotFoundError(f"No {filetype} file found for the passed in analysis: {ana.dataset_path}")
    else:
        # Otherwise, return the primary analysis for the dataset
        ana = get_primary_analysis(dataset_id, is_spatial)
    return ana

def get_primary_analysis(dataset_id, is_spatial=False) -> "SpatialAnalysis | Analysis":
    """Return the primary analysis for a dataset."""
    ds = Dataset(id=dataset_id, has_h5ad=1)
    filetype = "h5"

    # Ensure the zarr file is retrieved instead of the h5
    if is_spatial:
        ds.dtype = "spatial"
        filetype = "zarr"
        ds.has_h5ad = 0  # does not affect get_file_path but sanity-checking

    dataset_file_path = ds.get_file_path()

    # Let's not fail if the file isn't there
    if not Path(dataset_file_path).exists():
        raise FileNotFoundError(f"No {filetype} file found for this dataset {dataset_file_path}")
    if is_spatial:
        ana = SpatialAnalysis(type="primary", dataset_id=dataset_id)
    else:
        ana = Analysis(type="primary", dataset_id=dataset_id)
    return ana

def normalize_analysis_input(analysis: dict | str | None) -> dict | None:
    """
    Ensures that the analysis parameter is a dictionary if it is not None.

    Parameters
    ----------
    analysis : dict, str, or None
        The analysis parameter which can be a dictionary, a string (analysis ID), or None.

    Returns
    -------
    dict or None
        The analysis as a dictionary if it was provided as an string, otherwise returns it unchanged.
        If analysis is None, returns None.
    """
    if analysis is None:
        return None
    if isinstance(analysis, dict):
        return analysis
    if isinstance(analysis, str):
        return {"id": analysis}
    raise ValueError("Analysis must be a dict, str, or None.")

def test_analysis_for_zarr(analysis: "Analysis | SpatialAnalysis") -> bool:
    """
    Tests if the provided analysis is backed by a Zarr file store, in which AnnData objects are gotten in a different way.

    Parameters
    ----------
    analysis : Analysis or SpatialAnalysis
        The analysis object to be tested.

    Returns
    -------
    bool
        True if the analysis is a SpatialAnalysis and its type is 'primary', indicating it uses a Zarr file.
        False otherwise.
    """
    return isinstance(analysis, SpatialAnalysis) and analysis.type == "primary"

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

    def __init__(
        self,
        id=None,
        dataset_id=None,
        user_id=None,
        session_id=None,
        label=None,
        type=None,
        vetting=None,
    ):
        self.id = id
        self.dataset_id = dataset_id
        self.label = label
        self.session_id = session_id
        self.user_id = user_id

        if self.dataset_id is None:
            raise Exception("ERROR: Analysis object must have a dataset_id set.")

        if self.id is None:
            self.id = self.dataset_id

        # types are 'primary', 'public', 'user_saved', 'user_unsaved'
        self.type = type if type is not None else "primary"

        # vettings are None, 'owner', 'gear', or 'community'
        self.vetting = vetting

        # if user ID wasn't set but the session was, do the lookup
        if user_id is None and session_id:
            self.user_id = get_user_id_from_session_id(session_id)

        if self.type not in ["primary", "public", "user_saved", "user_unsaved"]:
            raise Exception(f"ERROR: Invalid type '{self.type}' for Analysis instance")

        if self.vetting not in [None, "owner", "gear", "community"]:
            raise Exception(f"ERROR: Invalid vetting '{self.vetting}' for Analysis instance")

    def __repr__(self):
        pipeline_file = self.settings_path
        with open(pipeline_file) as f:
            json_data = json.loads(f.read())
            # change "user_session_id" to "analysis_session_id" for consistency
            if "user_session_id" in json_data:
                json_data["analysis_session_id"] = json_data.pop("user_session_id")
            return json.dumps(json_data, indent=4)

    def _serialize_json(self):
        # Called when json modules attempts to serialize
        return self.__dict__

    @property
    def base_path(self) -> Path:
        """
        Returns the base directory path for an analysis, based on its actual type.  This allows for the support
        of parallel types 'primary', 'public', 'user_saved', 'user_unsaved'
        """
        if self.type == "primary":
            return self.primary_path

        else:
            # all other types require analysis ID to be set
            if self.id is None:
                raise Exception(
                    "ERROR: base_path called on Analysis object with no id attribute set."
                )

            if self.dataset_id is None:
                raise Exception(
                    "ERROR: base_path called on Analysis object with no dataset_id attribute set."
                )

            if self.type == "public":
                # ./$dataset_id/$analysis_id
                return root_dir / "www" / "analyses" / "by_dataset" / self.dataset_id / self.id

            elif self.type == "user_saved":
                if self.user_id is None:
                    raise Exception(
                        "ERROR: base_path called on Analysis object with no user_id attribute set. Probably not logged in."
                    )

                # ./$user_id/$dataset_id/$analysis_id/$dataset_id.h5ad
                return root_dir / "www" / "analyses" / "by_user" / str(self.user_id) / self.dataset_id / self.id

            elif self.type == "user_unsaved":
                if self.session_id is None:
                    raise Exception(
                        "ERROR: base_path called on Analysis object with no session_id attribute set. Probably not logged in."
                    )

                # /tmp/$session/$dataset_id/$analysis_id/$dataset_id.h5ad
                return Path(f"/tmp/{self.session_id}/{self.dataset_id}/{self.id}")

        raise Exception(f"ERROR: Invalid type '{self.type}' for Analysis instance")

    @property
    def dataset_path(self) -> Path:
        return self.base_path / f"{self.dataset_id}.h5ad"

    def discover_vetting(self, current_user_id: int | None = None):
        """
        This describes the public attribution of the analysis.  Making it a derived value via this
        method rather an than explicitly stored one in the database or JSON layer so that changes
        made to the dataset ownership won't require updating of this as well.  This method will
        just return the correct value.

        Returns one of the values 'owner', 'gear' or 'community'

        If the owner is also a curator it gives priority to the curator status (gear)
        """

        if self.type == "primary":
            # ? is this right
            self.vetting = "gear"

        # Analysis.user_id must be knownor we can't do this
        if self.user_id is None:
            raise Exception(
                "ERROR: Attempted to call Analysis.discover_vetting() without an owner assigned to the analysis"
            )

        current_user = get_user_by_id(current_user_id)

        if current_user is None:
            raise Exception(
                "ERROR: Attempted to call Analysis.discover_vetting() without a current user assigned to the analysis"
            )

        # if the user who created it is a gear curator, return that
        if current_user.is_curator:
            self.vetting = "gear"
        # Are the current user and dataset owner the same?
        elif self.user_id == current_user_id:
            self.vetting = "owner"
        else:
            self.vetting = "community"

        return self.vetting

    def discover_type(self)-> str | None:
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
        self.type = None

        type_checks = [
            ("primary", self.id == self.dataset_id),
            ("user_saved", (root_dir / "www" / "analyses" / "by_user" / str(self.user_id) / str(self.dataset_id) / str(self.id) / f"{self.dataset_id}.h5ad").exists()),
            ("user_unsaved", (Path(f"/tmp/{self.session_id}/{self.dataset_id}/{self.id}/{self.dataset_id}.h5ad").exists())),
            ("public", (root_dir / "www" / "analyses" / "by_dataset" / str(self.dataset_id) / str(self.id) / f"{self.dataset_id}.h5ad").exists()),
        ]
        for t, cond in type_checks:
            if cond:
                self.type = t
                break

        if self.type is None:
            raise Exception(
                f"ERROR: Unable to determine type for analysis {self.id} with dataset ID {self.dataset_id}."
            )
        return self.type

    @classmethod
    def from_json(cls, jsn):
        """
        Returns an Analyis object from a JSON object with the same attributes
        """

        # Dataset ID is now saved only in the dataset object, but some analyses may have it as a top-level attribute
        try:
            dataset_id = jsn["dataset"]["id"]
        except KeyError:
            # This is for legacy pipelines.  Both entries should be present though
            dataset_id = jsn["dataset_id"]

        try:
            session_id = jsn["analysis_session_id"]
        except KeyError:
            # This is for legacy pipelines.  This is the old name before it was realized
            # that a user could have analyses across multiple user sessions
            session_id = jsn["user_session_id"]

        ana = Analysis(
            id=jsn["id"],
            dataset_id=dataset_id,
            label=jsn["label"],
            session_id=session_id,
            user_id=None,
            type=jsn["type"],
        )

        ## get the rest of the properties
        for k in jsn:
            if not hasattr(ana, k):
                # some were manually named, skip them
                if k not in ["user_session_id", "analysis_session_id"]:
                    setattr(ana, k, jsn[k])

        return ana

    def get_adata(self, **kwargs) -> "AnnData":
        """
        Returns the anndata object for the current analysis.
        """
        dataset_path = self.dataset_path
        if "dataset_path" in kwargs:
            dataset_path = kwargs["dataset_path"]
            kwargs.pop("dataset_path")

        return H5adAdapter(dataset_path).get_adata(**kwargs)

    @property
    def marker_gene_json_path(self):
        return f"{self.base_path}/{self.dataset_id}.marker_gene_table.json"

    def _parent_path_by_type(self, atype=None) -> Path:
        """
        Returns the base directory base path for an analysis, based on a hypothetical type.  This allows for the support
        of parallel types 'primary', 'public', 'user_saved', 'user_unsaved'.

        This is generally useful if you want to find where an analysis directory would be if it existed, such as
        when searching for lists of analyses.
        """

        if atype == "primary":
            return root_dir / "www" / "datasets"

        else:
            if self.dataset_id is None:
                raise Exception(
                    "ERROR: _parent_path_by_type() called on Analysis object with no dataset_id attribute set."
                )

            if atype == "public":
                return root_dir / "www" / "analyses" / "by_dataset" / self.dataset_id

            elif atype == "user_saved":
                if self.user_id is None:
                    raise Exception(
                        "ERROR: _parent_path_by_type() called on Analysis object with no user_id attribute set."
                    )

                return root_dir / "www" / "analyses" / "by_user" / str(self.user_id) / self.dataset_id

            elif atype == "user_unsaved":
                if self.session_id is None:
                    raise Exception(
                        "ERROR: _parent_path_by_type() called on Analysis object with no session_id attribute set."
                    )

                return Path("/tmp") / self.session_id / self.dataset_id

        raise Exception(f"ERROR: Invalid type '{atype}' for Analysis instance")

    @property
    def primary_path(self) -> Path:
        return root_dir / "www" / "datasets"

    @property
    def settings_path(self) -> Path:
        """
        Returns the file path to the pipeline settings JSON file for the current dataset.

        The path is constructed by joining the base path
        and the dataset ID, resulting in a file named '<dataset_id>.pipeline.json'.

        Returns:
            str: The full path to the pipeline settings JSON file.
        """
        return Path(self.base_path) / f"{self.dataset_id}.pipeline.json"


class SpatialAnalysis(Analysis):
    """
    Spatial-based analysis object.  Inherits from Analysis and adds spatial-specific methods.
    """

    def __init__(
        self,
        id=None,
        dataset_id=None,
        user_id=None,
        session_id=None,
        label=None,
        type=None,
        vetting=None,
    ):
        super().__init__(
            id=id,
            dataset_id=dataset_id,
            user_id=user_id,
            session_id=session_id,
            label=label,
            type=type,
            vetting=vetting,
        )

        self.set_adapter()

    @property
    def dataset_path(self) -> Path:
        if self.type == "primary":
            return self.primary_path / f"{self.dataset_id}.zarr"
        else:
            return self.base_path / f"{self.dataset_id}.h5ad"

    def determine_platform(self, sdata):
        try:
            platform = sdata.tables["table"].uns["platform"]
            return platform
        except KeyError:
            raise ValueError("No platform information found in the dataset")

    def discover_type(self) -> str | None:
        super().discover_type()
        # This is a good time to set the adapter too
        self.set_adapter()
        return self.type

    def get_sdata(self, **kwargs) -> "SpatialData":
        """
        Returns a SpatialData object by loading it from a .zarr file.

        This method is only supported when the adapter class is ZarrAdapter. If the adapter class is not ZarrAdapter, a ValueError is raised.

        Returns:
            SpatialData: The loaded spatial data object.

        Raises:
            ValueError: If the adapter class is not ZarrAdapter.
        """
        if self.adapter_cls is not ZarrAdapter:
            raise ValueError("get_sdata() is only supported for .zarr files")

        adapter = self.adapter_cls(self.dataset_path)
        return adapter.get_sdata()

    def get_adata(self, **kwargs) -> "AnnData" :
        """
        Retrieve an AnnData object using the configured adapter class.

        This method initializes the adapter class with the dataset path and calls its
        `get_adata` method, passing any additional keyword arguments.

        Raises:
            ValueError: If no adapter class is set for this analysis.

        Returns:
            AnnData: The annotated data object retrieved by the adapter.
        """

        if self.adapter_cls is None:
            raise ValueError("No adapter class set for this analysis")

        adapter = self.adapter_cls(self.dataset_path)
        return adapter.get_adata(**kwargs)

    @property
    def primary_path(self) -> Path:
        return root_dir / "www" / "datasets" / "spatial"

    def set_adapter(self) -> None:
        """
        Determines and returns the appropriate adapter class based on the object's type.

        Returns:
            typing.Type[H5adAdapter] | typing.Type[ZarrAdapter] | None:
                - ZarrAdapter if self.type is "primary"
                - H5adAdapter if self.type is not "primary" and not None
                - None if self.type is None
        """
        if self.type == "primary":
            self.adapter_cls = ZarrAdapter
        elif self.type is None:
            self.adapter_cls = None
        else:
            self.adapter_cls = H5adAdapter

class AnalysisCollection:
    def __init__(self, public=None, user_saved=None, user_unsaved=None):
        self.public = [] if public is None else public
        self.user_saved = [] if user_saved is None else user_saved
        self.user_unsaved = [] if user_unsaved is None else user_unsaved

    def __repr__(self):
        return json.dumps(self.__dict__)

    def _collect_analysis_json(self, ana: "Analysis | SpatialAnalysis", atype: str) -> list:
        """
        Searches a directory to find stored analysis files.  Assumes the dir passed is a parent directory which
        can contain more than one analysis directory.  Looking essentially for this:

        $dir/*/*.pipeline.json

        Returns a list of JSON objects, one for each analysis
        """
        analyses = list()
        a_dir = ana._parent_path_by_type(atype=atype)

        ana_cls = type(ana)

        if a_dir.exists():
            if atype == "primary":
                ana.type = "primary"
                json_path = ana.settings_path

                if json_path.exists():
                    json_obj = json.loads(open(json_path).read())
                    analyses.append(ana_cls.from_json(json_obj))
            else:
                for json_thing in a_dir.glob("**/*.pipeline.json"):
                    json_path = json_thing
                    json_obj = json.loads(open(json_path, encoding="utf-8").read())
                    analyses.append(ana_cls.from_json(json_obj))

        return analyses

    @property
    def _serialize_json(self):
        # Called when json modules attempts to serialize
        return self.__dict__

    def get_all_by_dataset_id(self, user_id=None, session_id=None, dataset_id=None, is_spatial=False):
        """
        Gets all possible analyses for a given dataset, including primary, public, user-saved and
        user-unsaved analyses.
        """
        # clear any existing ones first
        self.__init__()

        ## Create hypothetical analysis to get paths
        ana_cls = Analysis
        if is_spatial:
            ana_cls = SpatialAnalysis
        ana = ana_cls(dataset_id=dataset_id, user_id=user_id, session_id=session_id)

        ## Each of these is a list of JSON objects
        self.primary = self._collect_analysis_json(ana, "primary")
        self.public = self._collect_analysis_json(ana, "public")
        if user_id:
            self.user_saved = self._collect_analysis_json(ana, "user_saved")
        if session_id:
            self.user_unsaved = self._collect_analysis_json(ana, "user_unsaved")

### File-specific adapters to load data from different formats ###

class H5adAdapter:
    """
    Adapter to provide a common interface for accessing AnnData objects from different sources.
    """

    def __init__(self, h5ad_path: Path):
        self.h5ad_path = h5ad_path

    def get_adata(self, backed=False, force_sparse=False) -> "AnnData":
        """
        Load and return the AnnData object from the specified H5AD file.

        Parameters:
            backed (bool, optional): If True, open the H5AD file in backed mode (read-only),
                which allows working with data that does not fit into memory. Defaults to False.
            force_sparse (bool, optional): If True, load the data matrix as a sparse array.
                Defaults to False.

        Returns:
            AnnData: The loaded AnnData object.

        Notes:
            - Requires the Scanpy library (`sc`).
            - The `h5ad_path` attribute must be set to the path of the H5AD file.
        """
        kwargs = {}

        import scanpy as sc

        if backed:
            kwargs["backed"] = "r"
        if force_sparse:
            kwargs["as_sparse"] = "raw.X"

        return sc.read_h5ad(self.h5ad_path, **kwargs)

class ZarrAdapter:
    """
    Adapter to provide a common interface for accessing SpatialData and AnnData objects from different sources.
    """

    def __init__(self, zarr_path: Path):
        self.zarr_path = zarr_path

    def get_sdata(self) -> "SpatialData":
        """
        Retrieves the SpatialData object from the specified Zarr path.

        Returns:
            SpatialData: The loaded SpatialData object from the Zarr file.

        Raises:
            FileNotFoundError: If the Zarr file does not exist at the specified path.
        """
        import spatialdata as sd

        if not self.zarr_path.exists():
            raise FileNotFoundError(f"Dataset not found at {self.zarr_path}")
        return sd.read_zarr(self.zarr_path)

    def get_adata(self, include_images=None) -> "AnnData" :
        """
        Generate and return an AnnData object from the current spatial dataset.

        This method retrieves the spatial dataset, determines its platform type,
        and uses the appropriate spatial handler to process the data. Optionally,
        images can be included in the resulting AnnData object.

        Args:
            include_images (bool, optional): Whether to include images in the AnnData object.
                If None, this is determined automatically based on the spatial handler.

        Returns:
            AnnData: The processed AnnData object, with additional metadata in `uns`:
                - "has_images": Whether images are available in the dataset.
                - "img_name": The name of the image, if available.

        Raises:
            ValueError: If platform information is missing or unsupported.
        """
        sdata = self.get_sdata()

        from gear import spatialhandler

        # Ensure the spatial data type is supported
        platform = None
        try:
            platform = sdata.tables["table"].uns["platform"]
        except KeyError:
            raise ValueError("No platform information found in the dataset")

        if not platform or platform not in spatialhandler.SPATIALTYPE2CLASS:
            raise ValueError(f"Invalid or unsupported spatial data type {platform}")

        spatial_obj = spatialhandler.SPATIALTYPE2CLASS[platform]()
        spatial_obj.sdata = sdata

        # Filter by bounding box (mostly for images)
        spatial_obj.filter_sdata_by_coords()

        if include_images is None:
            include_images = spatial_obj.has_images

        # Create AnnData object
        # Do not include images in the adata object (to make it lighter)
        spatial_obj.convert_sdata_to_adata(include_images=include_images)
        adata = spatial_obj.adata
        # Extra metadata to help with determining if images are available and where they are
        adata.uns["has_images"] = spatial_obj.has_images
        if spatial_obj.has_images:
            adata.uns["img_name"] = spatial_obj.img_name
        return adata
