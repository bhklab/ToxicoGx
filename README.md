# ToxicoGx
A R package in-progress. Goal is to analyze the toxicity of drugs across large-scale toxicogenomic datasets.

Order to run files:
1. curationDrugObject.R 
2. featureInfo.R
3. phenoInfo.R
4. esetFINAL.R 
5. cellObject.R 
6. drugObject.R 
7. curationCellObject.R 
8. curationTissueObject.R 
9. sensitivityObjects_DNA.R
10. sensitivityObjects_LDH.R
11. createTSet.R

- 1, 2 can be run interchangeably
- 3 must be run after curationDrugObject.R
- 4 must be run after esetFINAL.R
- 5-10 can be run interchangeably, but must be run after phenoInfo.R
- Choose which type of viability data is desired for the tSet (DNA(%) (9), LDH(%) (10))
- 11 must be run after all other files have been run
