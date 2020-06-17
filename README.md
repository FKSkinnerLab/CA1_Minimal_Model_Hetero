# CA1_Minimal_Model_Hetero

Anton Lunyov (uploaded June 2020), code from summer 2017+;
FKSkinner (edits)

Language: python

## Summary ##

This repository includes python code and scripts in 2 folders (1. Cell_Database_Metrics, and 2. Network_Sims) as used in the Chatzikalymniou et al. (2020) paper.  The code does the following: 

  1. Generates the model database (a,b,d,klow parameters), along with the generated model database output; obtains SFA, Rheo, PIR metrics from generated model database; creates histograms (as in paper Figure). 
  FILES are bubble... npy ETC.  The pyramidal cell models are used to create heterogenous populations to use in the E-I network simulations.
  
  2. Simulates heterogeneous E-I networks using brian2, with pyramidal cell models from the generated model database (generated npy files needed); analyzes output for power (FFT), as well as other aspects; plotting of output for voltages, IPSCs, EPSCs for some 
  FILES are BLAH.
  
  ## 3. Simulated network ouput of the heterogeneous output, and homogeneous ones. etc. (could go on OSF instead maybe? Anton?) ##

Izhikevich cellular model equations and parameters use the strongly adapting CA1 pyramidal cell model; E-I network simulations done via brian2 software (briansimulator.org) and network rationale and details can be found in Ferguson et al. 2017.  
Parameter values start from the 8th row of Table 5.  All details can be found in Chatzikalymniou et al.




## References ##

Chatzikalymniou, A. P.; Gumus, M.; Lunyov, A; Rich, S.; Lefebvre, J.; Skinner, F.K (submitted, 2020). 
Translating mechanisms from minimal to detailed models of CA1 hippocampal microcircuits reveals how theta rhythms emerge and how their frequencies are controlled - bioarxiv 

Ferguson, K.A., Chatzikalymniou, A.P., Skinner, F.K., 2017. Combining Theory, Model, and Experiment to Explain How Intrinsic Theta Rhythms Are Generated in an In Vitro Whole Hippocampus Preparation without Oscillatory Inputs. eNeuro 4. https://doi.org/10.1523/ENEURO.0131-17.2017

