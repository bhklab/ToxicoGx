# ToxicoGx
A R package in-progress. Goal is to analyze the toxicity of drugs across large-scale pharmacogenomic datasets.

Order to run files:
1. featureInfo.R
2. phenoInfo.R
3. esetFINAL.R 
4. cellObject.R 
5. drugObject.R 
6. curationCellObject.R 
7. curationDrugObject.R 
8. curationTissueObject.R 
9. sensitivityObjects_DNA.R
10. sensitivityObjects_LDH.R
11. createTSet.R

- 1, 2 can be run interchangeably
- 3-10 can be run interchangeably, but must be run after phenoInfo.R
- Choose which type of viability data is desired for the tSet (DNA(%) (9), LDH(%) (10))
- 11 must be run after all other files have been run
