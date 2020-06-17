# CA1_Minimal_Model_Hetero

Anton Lunyov (uploaded June 2020), code from summer 2017+;
FKSkinner (edits)

Language: python

## Summary ##

This repository includes python code and scripts that 

  i. generates the model database (a,b,d,klow parameters), along with the generated model database output (TO ADD FILE NAMES)
  
  ii. obtains SFA, Rheo, PIR metrics from model database

Izhikevich cellular model equations based on strongly adapting CA1 pyramidal cell model; E-I network simulations done via brian2 software.



Contains code for network simulations using the CellMetrics generated parameters

The Default network simulation contains the code used to run a heterogeneous network of neurons based on scenario 8 in Brian2, and the other versions of the code are based on it. The cell metrics generated for the runs are found in folders where they are needed to run as .npy files



## References ##

Chatzikalymniou, A. P.; Gumus, M.; Lunyov, A; Rich, S.; Lefebvre, J.; Skinner, F.K (submitted, 2020). 
Translating mechanisms from minimal to detailed models of CA1 hippocampal microcircuits reveals how theta rhythms emerge and how their frequencies are controlled - bioarxiv 

Ferguson, K.A., Chatzikalymniou, A.P., Skinner, F.K., 2017. Combining Theory, Model, and Experiment to Explain How Intrinsic Theta Rhythms Are Generated in an In Vitro Whole Hippocampus Preparation without Oscillatory Inputs. eNeuro 4. https://doi.org/10.1523/ENEURO.0131-17.2017

