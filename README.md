geoRcb
======

An extension of the `R`package [`geoR`](http://www.leg.ufpr.br/geoR/) that works with **cost-based distances**

Basically, three functions are adapted in order to admit custom distance matrices as additional arguments. 
Namely `variog`, to produce empirical variograms; `likfit`, to fit theoretical variograms and `krige.conv`, to preform kriging.
Check the documentation on these functions for examples of usage.

The cost-based distance matrix need to be computed elsewhere.
I have used [`GRASS GIS`](https://en.wikipedia.org/wiki/GRASS_GIS), with the help of a [script](https://github.com/famuvie/v.costdist.mat).
However, currently there is an R-package [gdistance](http://cran.r-project.org/web/packages/gdistance/index.html) which could help in computing cost-based distances.

## Instalation
  ```R
  install.packages('devtools')  # package devtools needed
  devtools::install_github('famuvie/geoRcb')
  ```

## Examples
   ```R
   library(geoRcb)
   example('variog')
   example('likfit')
   example('krige.conv')
   ```

## References 

**Antonio López-Quílez and Facundo Muñoz (2009)**. Geostatistical computing of acoustic maps in the presence of barriers.
*Mathematical and Computer Modelling* **50**(5-6):929-938.
[preprint](http://www.uv.es/famarmu/doc/preprint-mcm2009.pdf) | [journal](http://dx.doi.org/10.1016/j.mcm.2009.05.021)
