readme.txt
this file


normalizedLabelingData_fastedCD.xlsx
example input file containing the normalized labeling data (tissue data is not required if only calcualting direct contributions and fluxes between nutrients)


FcircPrime_fastedCD.xlsx
example input file containing the carbon atom Fcirc values (and errors) for the list of nutrients


calDirectTCAContributions.m
script to calcuate the direct TCA contributions from nutrients


calDirectContributions.m
script to calcualte the direct contributions between nutrients


calInterconvertingFluxes.m
script to calculate the absolute interconverting fluxes between nutrients. the input files are the direct contributions and the carbon atom Fcirc values


myvariance.m
a function used by the scripts "calDirectContributions.m" and "calDirectTCAContributions.m"
