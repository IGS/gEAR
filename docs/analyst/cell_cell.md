# Cell-cell communication

The cell-cell communication page allows for visualization of the output of a [dominoSignal](https://github.com/FertigLab/dominoSignal) analysis run on a chosen data.

Once the output of such a run is uploaded to gEAR (in the form of a .parquet file and an accompanying JSON object containing the relevant metadata), the cell-cell communication page can be used to interactively visualize the output.

Users can generate dot plots showing the ligand and receptor expression across chosen ligand-receptor and sending-receiving cluster pairs, with filters allowing them to narrow down to only interactions between selected groups or to interactions whose expression exceeds user-specified thresholds. The page also creates a series of river plots (also referred to as Sankey diagrams), showing the entire pathway from sending to receiving clusters, via the ligands, receptors and transcription factors included in the input DominoSignal object. As with the dot plots, users can filter these plots down to include only selected steps of the pathway, include only interactions that involve selected groups, or apply the same expression thresholds.

There are further options to adjust the sizes of elements of the plots, alter the background color, or adjust the download options for the generated figures.