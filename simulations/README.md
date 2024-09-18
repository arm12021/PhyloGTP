# Simulation Scripts Instructions

The scripts located in this directory can be used to recreate the simulated datasets used in the evaluations. First, the *genTrue.sh* script creates the true species tree and error-free gene trees. Second, the *genInferred.sh* script generates sequences for each of the gene trees generated in the previous step, then applies iqtree2 to compute ML error-prone gene trees. We also provided a MrBayes nexus file which is intended to create the gene tree samples from the sequences used to evaluate AleRax. Note that the sequences must first be converted from Phylip to nexus format, and the user should specify a starting tree in the input nexus file named 'usertree'. The script can then be executed with
```
mb -i genAleRaxInput.nex
```
Note that the leaf labels will needed to be specified in the input mapping file of AleRax.