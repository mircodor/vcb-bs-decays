# Vcb fitter for Bs -> Ds(*) mu nu decays
This repository contains the code that will allow the users to perform a fit to the data published by the LHCb collaboration in Phys. Rev. D101 (2020) 072004. This fit will additionally employ several theory and experimental inputs, such as:
* Bs -> Ds mu nu (f+ function)
  * HPQCD, PRD101 (2020) 7, 074513
  * MILC, PRD85 (2012) 114502
  * LCSR, EPJC80 (2020) 4, 347
* Bs -> Ds* mu nu (V, A1 and A2 functions)
  * HPQCD, arXiv:2105.11433
  * HPQCD, PRD99 (2019) 114512
  * LCSR, EPJC80 (2020) 4, 347
  * LHCb, JHEP12 (2020) 144

The fit can be performed using the CLN parameterisation or the BGL parameterisation (up to order 10 of the series).

## Running the code
The code has been test with `ROOT 6.22.06`. In order to run the code this version of `ROOT` must be used (earlier versions might work too).
The first thing to be done is to generate the signal templates that will be used inside the fit. This is done issuing the following command: `root -l runGeneration.C`. This command executes a function (`generateSignalTemplates`) that generates the signal templates. This function takes four arguments that can be customised:
* `filename`: a string specifying the configuration file that needs to be read by the script. This configuration file contains the form-factor parameters used during the generation of the signal templates;
* `nentries`: the number of events that will be generated for Ds and Ds* decays.
* `rateDS_4D`: set this to `true` if you want to generate the signal template for Ds* using the 4-dimensional decay rate (q2 and three helicity angles). If you set this value to `false`, the 1-dimensional decay rate will be used instead. **N.B:** the 4-dimensional generation may take a long time.
* `_res2D`: set this to `true` to use the 2-dimensional resolution for the p_perp observable available in the supplementary material of the paper at this [link](http://cds.cern.ch/record/2706102/files/). If set to `false`, the 1-dimensional resolution is used instead.

The acceptance due to the LHCb detector and selection will be automatically taken into account and applied after the data have been generated.

A file called `file_templates.root` that contains the signal templates for Ds and Ds* with 500k events is already added to the repository.
After the signal templates have been generated, one can lauch the fit via the command `root -l RunFit.C`. This macro has just one argument, `configFile`, that specifies which configuration file needs to be used for the fit.

## Configuration files
The configuration files are contained in the directory `config/` and have the following structure (taking as example the CLN parameterisation):

```
Templates file_templates.root templTree 
==========================
FF-Model-Ref-P  CLN 2 
G0     1.074
rho2P  1.27 
-----------------------
FF-Model-Ref-V  CLN 4 
F1     0.902
rho2V  1.23 
R1     1.34 
R2     0.83
==========================
FF-Model-Fit-P  CLN 2
G0     1.07   0.04  0 0  0
rho2P  1.23   0.05  0 0  0
--------------------------  
FF-Model-Fit-V  CLN 4
F1     0.902  0.013 0 0  1
rho2V  1.23   0.1   0 0  0  
R1     1.34   0.1   0 0  0 
R2     0.83   0.1   0 0  0
---------------------------------
Other-Parameters 15
Vcb        40e-3    0.1      0 0  0
etaEW      1.0066   0.0050   0 0  1
tauBs      1.510    0.       0 0  0
fsfdBrDs   0.01265  0.00053  0 0  1  
BrD        0.00993  0.00024  0 0  1
BrDst      0.323    0.006    0 0  1
BrBdD      0.0231   0.001    0 0  1
BrBdDst    0.0505   0.0014   0 0  1
NrefD      36.4e3   1.6e3    0 0  1
NrefDst    27.8e3   1.2e3    0 0  1  
effRD      1.568    0.008    0 0  1
effRDst    1.464    0.007    0 0  1
physBkg    15560    0.       0 0  0
combBkg     5940    0.       0 0  0
normLHCb   0.6812   0.1      0 0  0
---------------------------------
Correlations   1
NrefD NrefDst -0.703
---------------------------------
Theory-Inputs 5 HPQCD LCSRDsS LCSRDs LHCb-PAPER-2019-046 MILC

```
The first line (`Templates`) is formed by the name of the file containing the signal templates and by the name of the `TTree` that needs to be read. It should not be necessary to change it. The lines containing `FF-Model-Ref-V` and `FF-Model-Ref-P` specify the model used during the generation of the signal templates (CLN in both cases) and the number of parameters of the model. These parameters are set to the values measured by the LHCb collaboration in [Phys. Rev. D101 (2020) 072004](https://inspirehep.net/literature?sort=mostrecent&size=25&page=1&q=find%20eprint%202001.03225). 


The lines `FF-Model-Fit-P` and `FF-Model-Fit-V` indicate the model used during the fit to the data (CLN or BGL) and the number of parameters. The parameters that will be used are indicated via the following syntax `parName     parValue  parError parLowLimit parUpperLimit  isGaussianConstrained`, which specifies the parameter name, central value, uncertainty, lower and upper limits and if the parameter is Gaussian constrained in the fit (`1` if yes, `0` if not). 


The line `Other-Parameters` indicates the number of other parameters used in the fit.

The line `Correlations` is followed by the a number that indicates the number of correlated parameters used in the fit (excluding the form-factor parameters). The correlations are simply indicated using the name of both parameters followed by the value of the correlation coefficient.

The line `Theory-Inputs` contains the number and the names of the theory/experimental inputs that will be used in the fit. Available options are: `HPQCD_Ds`, `HPQCD_DsS`, `LCSRDs`, `LCSRDsS`, `LHCb-PAPER-2019-046` and `MILC`. The fit will automatically include these additional inputs and minimise the `chi2` simultaneously.




