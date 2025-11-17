
# Welcome to the gEAR documentation

In addition to the text based documentation here, additional information in other formats (video, slides) is available at [https://umgear.org/manual.html](https://umgear.org/manual.html).

## What is gEAR?

The gene Expression Analysis Resource (gEAR) is a website for visualization and analysis of multiomic data both in public and private domains. The gEAR enables upload, visualization and analysis of bulk RNA sequencing, scRNA-seq data, and epigenome data.

gEAR is under active development and we are continuously adding new features. If you need assistance with any issues or have suggestions, please contact our user helpdesk [https://umgear.org/contact.html](https://umgear.org/contact.html).

The gEAR is maintained by a team at the University of Maryland School of Medicine, and led by Dr. Ronna Hertzano and Joshua Orvis. The gEAR is supported by the NIDCD/NIH R01DC013817, R01DC019370, NIMH/NIH R24MH114815 and the Hearing Health Foundation (Hearing Restoration Project).

## gEAR 2.0 Changes

- New thematic redesign
- New navigation panel on the left hand sidebar (follows from page to page)
- Renaming of various items/pages
  - Single-gene and Multi-gene curators are now called Single-Gene and Multi-gene Displays
  - profiles are now called Dataset Collections
  - gene carts are now called Gene Lists
- Home page has been changed to "My Workspace" Dashboard
- Resume where you left off history when logged in from the My Workspace Dashboard
- More customizations on panel size for Dataset Collections
- Video guides on how to perform various actions are found on specific pages
  - (Coming Soon) Interactive step-by-step guides
- Citation is now included in the sidebar to easily copy to manuscripts

## Basic Features

### Creating a gEAR account

To create an account, click the "Log In" dropdown button at the top right of the homepage and click on "Sign Up" or follow the link here ([Create gEAR Account](https://umgear.org/create_account.html)). On this page you will enter some basic information about you and your email address (which can be used to recover your password).  You will be required to authenticate the email address you provide. After creating an account, you will have access to all the tools available on the gEAR. To see your user information or edit any information, click on your username in the top right. At any point you can select to have the plots configured for a colorblind mode.

### Home page

The gEAR homepage ([umgear.org](https://umgear.org/)) has changed to a "My Workspace" Dashboard. Links to the most commonly used tools and utilities on the gEAR. At any time on the site, you can return to the homepage by clicking the gEAR icon in the top left corner.

<img width="700" alt="Screen Shot 2022-04-25 at 8 35 42 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/ScreenshotDashboard.png">

### Navigation Panel

There is now a navigation panel on the left hand side of the gEAR.  This navigation panel will follow between pages and allow the user to quickly access different modules within the gEAR.  As a user becomes more familiar with the symbols on the navigation panel, they can choose to minimize the panel to only display the symbols (without any text) by clicking on the arrow in the upper right hand corner of the navigation panel.

<img width="200" alt="Screen Shot 2022-04-25 at 8 35 42 AM"
src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/ScreenshotNavBar.png">

### Searching for genes of interest

Searching for genes of interest in public datasets is one of the most common uses of the gEAR platform. A good place to start searching is on the homepage. To search for your gene of interest, put your gene name in the gene expression search bar that appears about mid-page on the My Workspace Dashboard and press search.

<img width="700" alt="Screen Shot 2022-04-25 at 7 48 24 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/ScreenshotGeneExpressionSearch.png">

For more information about searching for genes of interest and other basic functions in gEAR see:

- [Walkthrough Video](https://www.youtube.com/watch?v=sr_kvm7W4OE)


### Finding datasets

If you wish to look at results for a particular dataset or explore what datasets are available, you can either search by clicking on the "Explore Datasets" tab on the My Workspace Dashboard or on "Datasets" tab on the Naviagation panel.

Once in the dataset explorer you can refine your search by type of data, date, organism, etc. by choosing options on the left side of the page. To search for genes within a particular dataset, click on the list view button at the top of the page (under "View:" the three lines button), click on the eye symbol and enter your gene/genes of interest in the search box.

<img width="700" alt="Screen Shot 2022-04-25 at 9 40 16 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/ScreenshotDataSetExplorer.png">


### Gene annotations

Gene annotations can be viewed on the gEAR platform in several ways. After searching for a gene, you will see several annotation results at the top of the screen.

- **Annotation by organism:**
  - There will be a dropdown box to select which organism you wish to choose for annotation. In the example below (searching for Atoh1) currently shows the mouse annotation as highlighted. Choosing other organisms will provide the name of the gene in that organism (e.g. ATOH1 in humans) and update the external links to match the gene name listed after "Showing results for gene:".
- **External resource links:**
  - These options open a new window with the pages for the searched gene on other web resources. Note: If your search includes more than one gene, the links will be displayed for the gene highlighted under search results in your search.

- **Deafness gene annotation:**
  - On gEAR there is additional annotation information related to deafness (for common deafness related genes). Clicking on the icon will provide a list of phenotypes and links out to deafness specific resources.

- **Functional annotation:**
  - Clicking on the functional annotation text will expand the window to show additional information for the gene of choice including: GO terms, full names, and any aliases.

<img width="700" alt="Screen Shot 2022-04-25 at 8 58 21 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/ScreenshotAnnotation.png">

### Viewing dataset information

There are two ways of viewing the information associated with a dataset:
1. Within the dataset explorer, click the list view expand button (in the View: section the open box) located at the top of the page.

<img width="700" alt="Screen Shot 2022-04-25 at 12 33 48 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/ScreenshotDatasetInformation.png">

*Still needs to be implemented working on it -- below*
2. When viewing a dataset collection (e.g. after searching for a gene), open the dropdown menu and click info for a popup window with additional information to appear.

*Need to change screenshot*

<img width="700" alt="Screen Shot 2022-04-25 at 12 36 38 PM" src="https://user-images.githubusercontent.com/28451557/165134122-80917bbd-d0e1-499d-a6b8-f7fc2edcb3be.png">

### Dataset Collections

dataset collections are gEAR's way of collecting together datasets around a similar topic. The gEAR team maintains a number of curated dataset collections but any user can create their own dataset collection of datasets to explore.

#### Selecting a premade collection

To select a dataset collection you wish to search within, use the dropdown menu in the gene expression viewer to pick a dataset collection. You will see options for highlighted dataset collections, dataset collections you own, dataset collections you've recently visited, dataset collections from groups you are a part of, and dataset collections shared with you.

<img width="700" alt="Screen Shot 2022-04-25 at 12 46 25 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/ScreenshotDatasetCollection_FromGeneExpression.png">

The default dataset collection (Hearing) contains 10 datasets that offer a broad overview of gene expression while other dataset collections may be devoted to specific organisms (e.g. zebrafish) or collections of displays used in individual publications. dataset collections can be created and modified in the Dataset Explorer (Note: you can't directly modify dataset collections owned by others but feel free to make duplicates with any modifications you wish).

#### Creating a new collection

In the dataset explorer, open collection management (on the top portion of the left bar) and choose the plus icon to create a new dataset collection. Add a dataset collection name and click add. After you click add, you will see the new dataset collection name in the dropdown menu. To add datasets to this dataset collection, click on list view under View, click the plus sign shown in the entry for each dataset, to remove already added datasets use the minus icon. If you want to add the same dataset multiple times just click the plus icon (This is useful for displaying the same dataset in two different visuals).  To see what dataset collections are currently in the dataset collection, click the eye icon. If you wish to make another dataset collection your default on the homepage (i.e. a dataset collection other than the default Hearing dataset collection) click the blue thumbs up arrow.

<img width="700" alt="Screen Shot 2022-04-25 at 12 49 32 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/ScreenshotCreateNewCollection.png">

#### Arranging a collection

The layout of collections (e.g. how they are viewed when searching) can be altered in the dataset explorer.

First, select the collection you wish to alter using the Collection Management dropdown (note: you can only alter collections you own, but feel free to make duplicates of other dataset collections and alter however you wish).

Next click on the arrangement view button (top center of page looks like a series of dots). If you do not have permission to alter the layout of the dataset collection, this button will not be present.

By dragging, you can rearrange the order and size of each dataset within the dataset collection.  You can now move the panels to be anywhere on the page, and sizing can be to make it larger/smaller or wider/narrower. Make sure to hit save when finished.

<img width="700" alt="Screen Shot 2022-06-10 at 8 48 31 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/ArrangeDatasetCollection.png">

### Epigenome displays -- Under Construction

The gEAR team is working on updating this module.  As this module is completed we will add documentation to help users utilize the new page.

## Analysis tools

gEAR incorporates several tools for analyzing data with the goal of improving the reusability of genomic data. Graphical interfaces are built around common analysis tools (e.g. a Seurat pipeline for single-cell data) to make the tools more accessible for those without coding experience or the need to download data. The tools can be accessed from the navigation panel or by clicking on the analysis tools dropdowns in the dataset explorer found with each dataset.

<img width="200" alt="Screen Shot 2022-04-25 at 2 20 00 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/AnalysisToolsFromNavBar.png">
<img width="700" alt="Screen Shot 2022-04-25 at 2 20 00 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/AnalysisToolsFromDatasetExplorer.png">

### Comparison tool

The comparison tool is designed for comparing differences in gene expression (e.g. bulk RNAseq) between two groups of samples and testing for statistical differences in expression (Note: If the dataset does not have replicates, statistical tests will be unavailable. For more in depth instruction of how to use this tool, see the instructions on the tool page (Need to update link) ([Link](https://umgear.org/compare_datasets.html)).

To get started with the RNAseq comparison tool, select the dataset you wish to use in the dropdown menu and the condition you would like to compare on your x and y axis. Groups can be included or excluded from the analysis using the dropdown and toggle menus.  There is a final dropdown for selecting comparison options that is optional, but allows you to configure the comparison further. After clicking the plot button, you will notice that genes less than 2 fold change are omitted from the plot (to increase plotting speed), this setting can be changed by increasing or decreasing the fold change cutoff option on the previous page by clicking on "Edit Parameters". Similarly, there are options for colorizing or only including genes on the plot which reach a set threshold for significance between groups. In order to perform statistical testing, the dataset must have a column named "replicate" with the replicate ID number for each sample.

<img width="700" alt="Screen Shot 2022-04-26 at 2 41 17 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/CompareTool.png">

### Single Cell Workbench

The workbench for single cell RNAseq (scRNAseq) is designed to allow biologists meaningful access to single cell data, even with limited informatics training. The workbench begins by selecting a dataset for analysis, and then offers analysis tools following several standard pre-processing steps. It was designed to mimic the workflow of a standardized Seurat pipeline. For more in depth instruction of how to use this tool, see the instructions on the tool page ([Link](https://umgear.org/analyze_dataset.html)).

After opening the single cell workbench and choosing a dataset, there will be a choice whether to conduct a new analysis or start from a stored analysis. Stored analysis starts with precomputed results, either uploaded by the dataset owner or from an analysis already conducted and saved on the gEAR platform.

To start at the beginning of the single cell processing steps (i.e. without using a stored analysis), first pick initial filtering criteria to apply to the data. Next, choose 'Apply filter' before moving to the next step. To advance between steps, use the drop down menus for the specific step.  In many cases, the menus will auto-advance you.  In steps that do not auto-advance this is intentional as some steps may have plots you want to look further into prior to advancing. You can also see what step you are on using the top timeline-like bar.  Clicking on any green step will take you back to that step; however, there are various steps you can't move back to (filtering steps) as the dataset has been altered. In each step of the pipeline, there are general recommendations for criteria to select but these can be changed to better fit your dataset or question.


<img width="700" alt="Screen Shot 2022-04-26 at 2 10 37 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/SCWorkbench.png">

### Multigene curator

The multigene curator ([Link](https://umgear.org/multigene_curator.html)) allows for the creation of displays and plots involving more than one gene (e.g. heatmaps, volcano plots, quadrant plots, etc).

#### Steps to create a multigene display

1. Select a dataset to plot
2. Choose whether to create a new plot or clone a previously curated plot
3. Choose a plot type and which saved analysis to utilize
4. Add gene names of interest.
Genes can either be added manually or via a gene list. If you start typing a gene name, all genes that match and are present in the dataset will be displayed in the dropdown menu.  Note: For Multi-gene curator you must select more than one gene.
5. Select what groups and plot attributes to include.
Many of these options can be changed on the next page as well.
6. Create plot.
After creating the plot, the display can either be saved as a display for the dataset used (Save New Display), or downloaded as an image. Note: if you make a multigene display the default view for a dataset, it will only be used as the default view if 'Multigene display' is toggled on during a search. Also during this step, you can create a new gene list from the genes selected to be included in the plot.

<img width="700" alt="Screen Shot 2022-04-26 at 10 11 12 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/MultiGeneDisplays.png">

### Gene lists

Similar to how dataset collections are used on gEAR to collect datasets, gene lists are used to make comparisons of multiple genes easier across multiple analyses. gene lists have various functions throughout the site but their most common use is to create user defined lists of genes that can be  added into multigene plots or other plotting functions. For example, if you find a set of genes in one analysis (or from a manuscript) and want to view how they differ between treatment groups in another dataset, you could either manually add the gene names each time you wish to use them or store them together in a gene list and be able to quickly access and share your gene lists with others.

gene lists can be created manually, via upload, or directly from analysis outputs:

- **Manually/ via upload:**

Open the gene list manager by clicking Gene Lists on the navigation panel on the left side of the page. The gene list Manager follows a similar layout as the dataset explorer where you can search for other public or private gene lists. To create a new list, click on the create new gene list button on the left side of the window. To upload from a file, choose the upload gene list option and select the file you wish to upload (file must be comma or space separated gene symbols). To add genes (either by copy/paste or manual entry, choose paste gene list. In addition to the genes you wish to add, give the list a name, description, and select the organism and privacy you wish for the cart (gene lists are private by default).

- **More information coming soon on weighted gene lists**

<img width="700" alt="Screen Shot 2022-04-25 at 3 01 02 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/GeneListManager.png">

<img width="700" alt="Screen Shot 2022-04-25 at 3 21 43 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/CreateNewGeneList.png">

- **From the RNAseq comparison tool:**

In the RNAseq comparison tool, you can create a new gene list containing your genes of interest directly from the plotting or statistical results. After making a plot, click and drag to select genes of interest.  Click on View and save selected genes to proceed. To download a text file listing all the genes use the download button (found next to save)) while the save button adds all the selected genes to a new gene list (You must provide a name to proceed) which can be used elsewhere.

<img width="700" alt="Screen Shot 2022-04-26 at 10 51 31 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/CompareToolSelectGenes.png">


<img width="500" alt="Screen Shot 2022-04-26 at 10 51 31 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/CompareToolMakeGeneList.png">

- **Using created or public gene lists**

gene lists can be used in the multigene display and gene expression tools. To use in the multigene viewer (and create your own visualiations), choose the dataset you wish to view and then under , use the dropdown menu "Quick search using Gene Lists" under select genes to select a gene list or individual genes from the lists to add your genes.

<img width="700" alt="Screen Shot 2022-04-25 at 3 30 54 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/UsingGeneListInMultiGeneCurator.png">

gene lists can be used on already created visualizations using the gene expression tool on the My Workspace Dashboard. First, click the Multi-gene Display radio button to change the search type from single to multigene. Then click on the "Quick search using Gene Lists" to select what gene list you would like to add. Press search to search for gene list visualizations of your genes of interest (Note: not all dataset collections have multigene visualizations, the multi-gene Displays tool can be used to create new visualizations).

<img width="700" alt="Screen Shot 2022-04-25 at 3 30 54 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/SearchDataSetCollectionForMultiGeneDisplay.png">

<img width="700" alt="Screen Shot 2022-04-25 at 3 30 54 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/SearchGeneExpressionUsingGeneList.png">

Multigene displays can also be toggled on the results page:

<img width="700" alt="Screen Shot 2022-04-25 at 3 42 29 PM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/ToggleToMultiGeneDisplaysWithGeneList.png">

## Upload data

For more information on how to upload different types of data into gEAR, see our [upload documentation](https://github.com/IGS/gEAR/blob/main/docs/Documentation/UploadingOverview.md).
Currently gEAR supports the upload of:

- Bulk RNAseq data
- Single Cell sequencing data
- Microarray data
- Epigenome data (Through integration with Gosling)

## Curate data/ build plots

The ability to create a variety of visualizations for the same dataset is one of the major strengths of the gEAR platform. In addition to the displays created by the initial uploader, if there are visualizations that are better suited to your question about a dataset you can quickly create new custom displays based on either your own analysis or primary analysis done by the uploader.

Besides for visualizations created in the Multi-Gene Displays curator or analysis tools (see above), custom displays can be added to dataset collections. To open a dataset in the curation tool (the same tool you see following initial dataset upload), click the Single-Gene or Multi-gene Displays on the left hand navigation bar and choose the dataset you wish to make a new display for and select create new display on the second step.  Continue on each step selecting options necessary for your visualization.

<img width="700" alt="Screen Shot 2022-04-26 at 10 19 28 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/OpenCuratorFromGeneExpression.png">

You can see displays shared and owned by you for a dataset when in a dataset collection.  First select the 3 dots at the upper right section of the dataset, and click on choose display. In the curation display window, you will see all the displays that are public or owned/shared with you. The current default view will be denoted greying out the selection "Set to Default". To curate a new display, click on curate new display at the top of the display selection dropdown and continue with other options.

<img width="700" alt="Screen Shot 2022-04-26 at 10 23 05 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/Single-GeneCurator-FromGeneExpression.png">

Choosing one of the curators will navigate to the curator tool. You will first have to choose whether to curate a new display or clone a previously curated display.

You will then have two choices for Dataset Type (Primary Data and Stored Analysis). Primary Data allows for plotting of data directly from the information uploaded (i.e. UMAP/tSNE coordinates or counts provided by the uploading individual). Choosing stored analysis allows for plotting from data that has been processed from the raw data in some way. Depending on the dataset, this could be an analysis uploaded by the dataset authors or a saved analysis you have conducted in one of the analysis tools (e.g. single cell workbench). Note: not all datasets will have a stored analysis (but you can create one).

After adding the analysis type and type of plot you wish to curate, you'll have other options based on your plot-type.  Fill these out as necessary and click plot where you can refine any options and see your visualization.

<img width="500" alt="Screen Shot 2022-04-26 at 10 31 33 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/SingleGeneCurator.png">

### Editing dataset information after upload

The dataset information shown in the dataset explore can be edited by the dataset's owner following upload. To edit, first click on the list view option at the top of the page.  Then click the pencil icon next to the dataset you wish to edit. On this page you can edit the dataset title, description, as well as add a Pubmed or GEO id. Note: Special characters cannot be used in the dataset title or descriptions.

<img width="700" alt="Screen Shot 2022-07-15 at 8 20 16 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/EditDataSetInfoAfterUpload.png">

## Downloading data

Public data in the gEAR platform can be downloaded in two ways:

- As a bundled file, containing the data files exactly as they were uploaded by the dataset owner. For more information on data formats used during upload, see our upload documentation ([Link](https://github.com/IGS/gEAR/blob/main/docs/Documentation/UploadingOverview.md))
- Via a formatted h5AD file, where each dataset upload is converted into an identical format. For more information on the h5AD file format see ([Link](https://anndata-tutorials.readthedocs.io/en/latest/getting-started.html))

To download a dataset, search a gene in the public dataset collection of interest then use the dropdown menu to select which download option you prefer.

<img width="700" alt="Screen Shot 2022-06-09 at 9 31 35 AM" src="https://github.com/IGS/gEAR/blob/documentation/docs/Documentation/Screenshots/DownloadData.png">
