# iPSCsInsulinResistance

## License
All code in this repository is subject to a GNU GPL v3. License. Please see LICENSE for the full license text.

## Description
This repository contains R code from a study investigating insulin resistance pathways in iPSCs.

Scripts files:
* residual expression computation [**RNAseqDataExplore.GENESIPS.IRvsIS.v6.317.88.short.R**]
* differential expression analysis [**RNAseqDataExplore.GENESIPS.IRvsIS.v6.317.88.short.R**, and for visualisation, **plotDEheatmaps.R**]
* co-expression analysis [**redo_coexp_IR-IS_v2-v3.R** which in turns calls **rnaSeqLimmaPlotPCAFunction_general.R** and **coexppMinervaSpecifiableParameters.R**]
* analysis of the atorvastatin experiment data [**GENESIPS_IRIS_RNAseq_Validation_combineData.R** for combining the featureCount processed data and **GENESIPS_IRIS_RNAseq_Validation_normaliseAdjustDE.R** for the analysis]
* the files **msigDB_enrichment_modified_subroutineOnly.R** and **usefulFunctionsForRNAseqAnalysis.R** are helper scripts called by some of the above.
