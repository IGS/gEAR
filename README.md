# gEAR Portal

**Description**:  The gEAR Portal was created as a data archive and viewer for gene expression data including microarrays, bulk RNA-Seq, single-cell RNA-Seq and more.  Initially created for the [hearing research community](https://umgear.org), instances of the gEAR have been cloned for other communities such as [brain research](https://nemoanalytics.org).

Other things to include:

  - **Technology stack**: The gEAR software is a LAMP stack utilizing Python and MySQL, H5AD for large expression matrix storage, D3 and Plot.ly for data visualization and has embedded [Epiviz](https://epiviz.github.io/) support to display epigenetic data.
  - **Status**:  This project has been in production for several years, though is in constant development so bugs certainly exist.
  - **Production / Demo instances**
	  - [UMgEAR](https://umgear.org) - Portal for hearing research
	  - [NeMO Analytics](nemoanalytics.org) - Portal for brain research

**Screenshot**: Example of home page after searching for a gene:

 ![](https://github.com/IGS/gEAR/blob/269f8f971301c15b69c50f3d11ad3441b2d24c78/docs/gear_overview.png)


## Installation

Setting up your own portal is admittedly a bit of work.  There are a lot of components to the portal and running on a server with at least 16 cores and 100GB+ of RAM is recommended.  The process is documented in the [setup.new_server.notes.md](docs/setup.new_server.notes.md) document.

## Usage

To learn how to use the software you can go to any existing portal and click the documentation link at the top, which will take you to a page like [this one](https://umgear.org/manual.html).  It has walk-through slides and YouTube videos for most topics.

## Known issues

This repository was hosted private for many years and we have just recently transitioned to this new public one.  If you find any issues, please use the Issue tracker here in GitHub to see if the problem is already reported or create your own ticket.

## Getting help

There are a few ways to get help with a gEAR Portal.  

 - Check the provided [documentation](https://umgear.org/manual.html).
 - Use the contact form at the top of any existing portal.  [Here](https://umgear.org/contact.html), for example.
 - Submit a ticket on the [issue tracker](https://github.com/IGS/gEAR/issues)


## Getting involved

If you'd like to contribute to the gEAR in any form (documentation, bug fixes, new features, etc.) just fork the project and create a pull request.


----

## Open source licensing info
The gEAR is distributed under the GNU AFFERO GENERAL PUBLIC LICENSE V3   For more info, see the [LICENSE](LICENSE) document.

[![DOI](https://zenodo.org/badge/289995740.svg)](https://zenodo.org/badge/latestdoi/289995740)



