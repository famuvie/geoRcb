## TODO prior to CRAN release:

There currently are a few violations to CRAN policy
http://cran.r-project.org/web/packages/policies.html

* Altering the namespace of another package should only be done with the
agreement of the maintainer of that package.

  We need to get the agreement of Paulo J. Ribeiro Jr <paulojus at ufpr.br> or
  (better) find a way to import the whole Namespace from geoR
  
* Packages should not modify the global environment (user’s workspace). 

  Instead of using a global flag .personal.definition.of.distances, we should
  find a more encapsulated way of flagging.

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
