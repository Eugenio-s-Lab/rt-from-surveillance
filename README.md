# rt-from-surveillance

This library implements the method to estimate the reproduction ratio of an infectious disease epidemic from surveillance data, accounting for the bias induced by the spatial structure of the epidemic. This method is described in the following publication:

Birello, P., Re Fiorentin, M., Wang, B., Colizza, V. & Valdano, E.
_Estimates of the reproduction ratio from epidemic surveillance may be biased in spatially structured populations._
Nature Physics 20, 1204-1210 (2024).
DOI: 10.1038/s41567-024-02471-7;
Open-access link: https://rdcu.be/dFJgn


Required inputs for the function 'corrected_rt_from_surveillance' are:
* a time series of reported or estimated infections in a number of spatial patches the geographic region under study is divided into;
* a co-location matrix, that accounts for the mobility patterns between the spatial patches;
* the population of spatial patches;
* generation interval distribution and distribution of delay from infection to reporting.

Further optional inputs are available.

To use our library, download the entire repository. Add your code to the repository, and *from surveillance_correction import corrected_rt_from_surveillance*. 
An example code and examples of input files can be found in this repository, in order to facilitate the use of the library.
