# CA1_Minimal_Model_Hetero

Anton Lunyov (uploaded June 2020), code from summer 2017+;
FKSkinner (edits)

Language: python

## Summary ##

This repository includes python code and scripts in 2 folders (1. Cell_Database_Metrics, and 2. Network_Sims) as implemented in the Chatzikalymniou et al. (2020) paper.  The code does the following: 

## Folder 1: Cell_Database_Metrics

Contains code that generates the model database (a,b,d,klow parameters), along with the generated model database output; stores the result in 3 4-D tensors corresponding to SFA, Rheo, PIR metrics from generated model database; plots histograms (as in paper Figure). The 4-D tensors are stored in files entitled adaptation_across_4.npy, rheo_across_4.npy and TC_across_4.npy. These files store the models that are subsequently used in the network simulations. 
  
  FILES 
  
  adaptation_across_4.npy - stores the SFA data as (initial slope - final slope) of the current-frequency curve across 4 varying dimensions (a,b,d,kLow)
  
  adaptation_across_4_ratio.npy - stores the SFA data as (initial slope / final slope) of the current-frequency curve across 4 varying dimensions (a,b,d,kLow)
  
  rheo_across_4.npy - stores the rheobase current data across 4 varying dimensions (a,b,d,kLow)
  
  TC_across_4.npy - stores the transition current data as a measure of post inhibitory rebound across 4 varying dimensions (a,b,d,kLow)
  
  gen_PIR_tensor.py - generates the post-inhibitory rebound data and stores it in TC_across_4.npy
  
  gen_rheo_tensor.py - generates the rheobase current data and stores it in rheo_across_4.npy
  
  gen_SFA_tensor.py - generates the spike frequency adaptation data and stores it in adaptation_across_4.npy and adaptation_across_4_ratio.npy
  
## Folder 2: Network_Sims 
  Contains code that simulates heterogeneous E-I networks adapted from Ferguson (2017) paper using brian2, with pyramidal cell models from the generated model database (the aforementioned .npy files); calculates and stores various output including power (FFT), plots of output for voltages, IPSCs, EPSCs, etc. 
  FILES 


Izhikevich cellular model equations and parameters use the strongly adapting CA1 pyramidal cell model; E-I network simulations done via brian2 software (briansimulator.org) and network rationale and details can be found in Ferguson et al. 2017.  
Parameter values start from the 8th row of Table 5.  All details can be found in Chatzikalymniou et al.



## References ##

Chatzikalymniou, A. P.; Gumus, M.; Lunyov, A; Rich, S.; Lefebvre, J.; Skinner, F.K (submitted, 2020). 
Translating mechanisms from minimal to detailed models of CA1 hippocampal microcircuits reveals how theta rhythms emerge and how their frequencies are controlled - bioarxiv 

Ferguson, K.A., Chatzikalymniou, A.P., Skinner, F.K., 2017. Combining Theory, Model, and Experiment to Explain How Intrinsic Theta Rhythms Are Generated in an In Vitro Whole Hippocampus Preparation without Oscillatory Inputs. eNeuro 4. https://doi.org/10.1523/ENEURO.0131-17.2017

