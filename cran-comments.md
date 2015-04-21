## Test environments
* Linux Mint 17 Qiana, R 3.2.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There were 3 NOTEs:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘"Facundo Muñoz" <facundo.munoz@orleans.inra.fr>’
New submission

  This is my first submission to CRAN

* checking foreign function calls ... NOTE
Foreign function calls to a different package:
  .C("binit", ..., PACKAGE = "geoR")
  .C("tgangle", ..., PACKAGE = "geoR")

  Actually, this package "implants" a modified version of a few functions
  from the dependency "geoR". Then everything happens within the Namespace
  of "geoR".
  
  This is done to avoid duplicating most of the code from geoR,
  while providing only the modified bits. The package geoR remains untouched,
  and it can always be loaded and used as expected. However, when geoRcb
  provides some additional functionality.
  
* checking R code for possible problems ... NOTE
Found the following possibly unsafe calls:
File ‘geoRcb/R/geoRcb.R’:
  unlockBinding(x, ns)
  unlockBinding(x, pe)

  This is the code needed to implant the extensions to the geoR package.
  

## Downstream dependencies
There are currently no downstream dependencies for this package
