# ProjectR scratch notes

## NeMO Archive

* Holds the weighted gene lists
  * Can grab some from other DB sources (believe this is lower priority)
* Holds tab_counts files (ROWmeta, COLmeta, DataMTX) which can be converted to H5AD
  * Should we still bundle and transfer to NeMO Analytics

## NeMO Analytics

Keep in mind this needs to ideally work on any gEAR portal

### API

* Run ProjectR
  * Loading input (all dataframe rows are genes)
    * Past analyses within datasets
      * Col - PCA, tSNE, UMAP
      * Stored in adata.obs if user uploaded or adata.obsm with "X_" prefix if performed in sc-RNAseq analysis workbench
      * ? Are there other keywords to look for?  Should these be selectable before calling the API (to display a specific combo)?
      * ? Does projectR even need to be run?  All stored PCA/tSNE/UMAP values are with respect to the observations, like in COLmeta_DIMRED files
    * Pattern repository
      * Col - PCs, Patterns, etc.
      * Currently I have patterns saved at /var/www/patterns
    * Weighted Gene carts
      * Col - Unnamed Weights
      * These are stored in the MySQL DB
  * Target dataset Input
    * Genes as rows
    * Observations (adata.obs) as cols
  * Projection patterns output dataframe
    * Rows are pattenrs
    * Cols are observations
    * Use adata.obs to tack on conditions, then can use those various conditions as params in the plot
* Create plots
  * ? What type of plots (scatter, heatmap, etc.)?
    * ? Single-pattern
    * ? Multi-pattern
  * ? Run through API calls or use lib/gear code?  Will I need to refactor?

### Gene aliases

* MySQL db table entry (or revival of previous one)
* Needed for cross-species projection
* Also useful for GCID work

### ? Projection curation page

* Currently building API code on projection.html but will this page be fleshed out for actual saved curations?
* This could get complicated if single- and multi-pattern curations have to be considered, especially if they have to go through the plotly_data and mg_plotly_data.py API calls
* Save these to database
  * Save pattern source (dataset/pattern/genecart)
  * ? Do we save the projection pattern dataframe in the target dataset's AnnData file for future quick-reference?
  * ? Can we toggle these configs back to "gene" mode (where genes are subbed for patterns in the config)

### Front page displays

* Toggle to switch between "genes" mode and "projection" mode
* ? I guess we add the "loadings" source and pattern selection options from projection.html into the sidepanel on this page?
* I'm guessing we would use the same default display configs for each dataset but substitute patterns for genes
  * This only applies if a default display exists, and the config is a plotly-based config (though I guess we can use Scanpy for tSNE/UMAP... just more complex)
* Separate displays for single-pattern and multi-pattern
* This UI may be in flux if the page gets redesigned

### Compare tool (lower priority)

* Need more infomation here

### Crons

* Perform Dimension Reduction under various analyses
  * PCA
  * CoGAPS
* ? Does this even belong here... should it go on NeMO Archive instead?