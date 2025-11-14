# Release test plan

This document describes the steps needed to perform a test on the back-end and interface
elements of gEAR before making a release into production.

## Server-side tests

None yet written. The /tests directory does have a collection of Excel sheets which can
be used to test the interface.

## Interface tests

- [Main page gene search](#main-page-gene-search)
- Profile creation, switch to primary
- Profile switch testing
- Dataset comparison tool
- sc-RNA analysis workbench
- Dataset upload
- Dataset curation
- Account creation

### Main page gene search

With the default profile loaded, search some canonical genes such as 'Pou4f3' and 'Sox2'
to make sure both display. Also test space and comma-separated strings.

