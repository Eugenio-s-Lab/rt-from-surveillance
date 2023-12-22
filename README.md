# rt-from-surveillance

This library is meant to implement the method to estimate the reproduction number of an infectious disease proposed in *Estimates of the reproduction ratio from epidemic surveillance may be biased in spatially structured populations*, Birello, Re Fiorentin, Wang, Colizza, Valdano. Preprint: https://arxiv.org/abs/2307.13798.

Required inputs for the function 'corrected_rt_from_surveillance' are:
* a time series of reported or estimated infections in a number of spatial patches the geographic region under study is divided into;
* a co-location matrix, that accounts for the mobility patterns between the spatial patches;
* the population of spatial patches;
* generation interval distribution and distribution of delay from infection to reporting.

Further optional inputs are available.

To use our library, download the entire repository. Place your code in the repository, and *from surveillance_correction import corrected_rt_from_surveillance*. 
An example code and examples of input files can be found in this repository, in order to facilitate the use of the library.
