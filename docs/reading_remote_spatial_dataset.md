# Reading remote spatial datasets

This is my attempts at reading in spatial datasets from a remote location.

## Apache config

One problem with reading Zarr from a remote area is that it is a directory store, not a file.  If you enable Indexing in the directory, it will show listings of files, but also create an HTML representation that has some prettied up elements.  If you curl the Zarr directory, you will get this HTML back.

When the `spatialdata.read_zarr` package is reading remote, it is going to rely on this HTML output for determining what groups are present.

This is the Apache entry that I got to work. I do not understand why ShowForbidden makes the listing HTML plain or if there are better alternatives to do the same thing but it works.

```
  <Directory /var/www/datasets/spatial/*.zarr>
    Options Indexes FollowSymLinks MultiViews
    Require all granted
    IndexOptions ShowForbidden
  </Directory>
```

## Spatialdata reading

### Method 1 - Reading as is

```python
import spatialdata as sd
rem_path = "https://devel.umgear.org/datasets/spatial/11692b64-b34a-4dbe-adc9-784a87a7a856.zarr"
sdata = sd.read_zarr(rem_path)
```

This results in the "shapes" group throwing an error as pyarrow cannot read the remote shapes.parquet file

### Method 2 - passing in selectors (best method, but still has some issues)

```python
import spatialdata as sd
rem_path = "https://devel.umgear.org/datasets/spatial/11692b64-b34a-4dbe-adc9-784a87a7a856.zarr"
sdata = sd.read_zarr(rem_path, selection=["images", "labels", "points", "tables"])
```

This works, but printing the REPR of sdata throws an error `AttributeError: 'NoneType' object has no attribute 'store'`

However you can access sdata properties like sdata.path, sdata.images, sdata.tables, the latter two are really what we need.

I tested this in the Panel dashboard, and the `spatialdata.to_legacy_anndata` function failed because it could not find the `sdata.tables["table"]` object.  It seems that if there are two tables in the .zarr (in this case "square_008um" and "table") one of them gets nested into the other when the spatialdata object is read.

We do not use "shapes" currently for our needs, but should we ever use "shapes", then we need to resolve the error from Method 1.

### Method 3 - specifying the store type

```python
import fsspec
import spatialdata as sd
rem_path = "https://devel.umgear.org/datasets/spatial/11692b64-b34a-4dbe-adc9-784a87a7a856.zarr"
store = fsspec.get_mapper(rem_path)
sdata = sd.read_zarr(store, selection=["images", "labels", "points", "tables"])
```

This leads to `TypeError: argument should be a str or an os.PathLike object where __fspath__ returns a str, not 'FSMap'` in trying to set the sdata.path property

Also tried the following:

```python
import zarr
root = zarr.open(rem_path)
sdata = sd.read_zarr(root, selection=["images", "tables"])
# or
sdata = sd.SpatialData.read(root, selection=["images", "tables"])
```

The TypeError from earlier returns but instead of FSMap it is a Group type
