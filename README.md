# ProxECAT
Proxy External Controls Association Test

rare-variant genetic region association test to robustly use publicly available data with allele frequencies in case-control analysis.

Notes:
--We have found better results by using fairly stringent internal minor allele frequency filters of 1% or 0.1% in either cases or in controls
--make sure to use cases and controls matched by ancestral population

To install in R, 
(1) install the package devtools
(2) library(devtools)
(3) install_github("hendriau/ProxECAT")



To run proxecat or proxecat_over in R,

proxecat::proxecat(x1, x2, x3, x4)
proxecat::proxecat_over(x1, x2, x3, x4, size=1000)

where 
x1	is a vector of the number of functional rare alleles in cases for each gene region

x2	is a vector of the number of synonymous rare alleles (or other proxy) in cases for each gene region

x3	is a vector of the number of functional rare alleles in controls for each gene region

x4	is a vector of the number of synonymous rare alleles (or other proxy) in controls for each gene region

size	is the overdispersion parameter. The smaller size, the more over dispersed. The larger size the less overdispersed. Size of >1000 is approximately Poisson.
