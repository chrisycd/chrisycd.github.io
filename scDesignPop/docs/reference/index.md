# Package index

## All functions

- [`calcParaVectors()`](https://chrisycd.github.io/scDesignPop/reference/calcParaVectors.md)
  : Generic function to compute model parameter vectors
- [`calcParaVectors(`*`<gam>`*`)`](https://chrisycd.github.io/scDesignPop/reference/calcParaVectors.gam.md)
  : A calcParaVectors Method for gam (mgcv package) Objects
- [`calcParaVectors(`*`<gamlss>`*`)`](https://chrisycd.github.io/scDesignPop/reference/calcParaVectors.gamlss.md)
  : A calcParaVectors Method for gamlss Objects
- [`calcParaVectors(`*`<glmmTMB>`*`)`](https://chrisycd.github.io/scDesignPop/reference/calcParaVectors.glmmTMB.md)
  : A calcParaVectors Method for glmmTMB Objects
- [`checkVectorContain()`](https://chrisycd.github.io/scDesignPop/reference/checkVectorContain.md)
  : Check membership of first vector compared to other vectors
- [`checkVectorEqual()`](https://chrisycd.github.io/scDesignPop/reference/checkVectorEqual.md)
  : Check if multiple vectors have same elements
- [`constructDataPop()`](https://chrisycd.github.io/scDesignPop/reference/constructDataPop.md)
  : Construct a list of input data
- [`constructDesignMatrix()`](https://chrisycd.github.io/scDesignPop/reference/constructDesignMatrix.md)
  : Construct a Design Matrix Dataframe
- [`constructEqtlGeno()`](https://chrisycd.github.io/scDesignPop/reference/constructEqtlGeno.md)
  : Construct eQTL annotation and genotype dataframe
- [`constructFormula()`](https://chrisycd.github.io/scDesignPop/reference/constructFormula.md)
  : Constructs a Model Formula
- [`constructPAFormula()`](https://chrisycd.github.io/scDesignPop/reference/constructPAFormula.md)
  : Construct a model formula for power analysis
- [`constructSCE()`](https://chrisycd.github.io/scDesignPop/reference/constructSCE.md)
  : Construct a SingleCellExperiment object with specified covariates.
- [`createPbulkExprGeno()`](https://chrisycd.github.io/scDesignPop/reference/createPbulkExprGeno.md)
  : Creates pseudobulk expression grouped by eQTL genotypes.
- [`example_eqtlgeno`](https://chrisycd.github.io/scDesignPop/reference/example_eqtlgeno.md)
  : Example genotype data for cell-type-specific cis-eQTLs
- [`example_eqtlgeno_Bcell`](https://chrisycd.github.io/scDesignPop/reference/example_eqtlgeno_Bcell.md)
  : Example B cell eQTL genotype data
- [`example_eqtlgeno_trans`](https://chrisycd.github.io/scDesignPop/reference/example_eqtlgeno_trans.md)
  : Example genotype data for cell-type-specific trans-eQTLs
- [`example_genopc_new`](https://chrisycd.github.io/scDesignPop/reference/example_genopc_new.md)
  : Genotype principal components for new individuals
- [`example_genopc_train`](https://chrisycd.github.io/scDesignPop/reference/example_genopc_train.md)
  : Genotype principal components for training individuals
- [`example_sce`](https://chrisycd.github.io/scDesignPop/reference/example_sce.md)
  : Example single-cell RNA-seq data (multi–cell type)
- [`example_sce_Bcell`](https://chrisycd.github.io/scDesignPop/reference/example_sce_Bcell.md)
  : Example B cell single-cell RNA-seq data with pseudotime
- [`extractFromSCE()`](https://chrisycd.github.io/scDesignPop/reference/extractFromSCE.md)
  : Extract data from SingleCellExperiment object
- [`extractFromSeurat()`](https://chrisycd.github.io/scDesignPop/reference/extractFromSeurat.md)
  : Extract data from Seurat object
- [`extractParaPop()`](https://chrisycd.github.io/scDesignPop/reference/extractParaPop.md)
  : Extract parameter matrix for a new covariate data frame
- [`fitCopulaPop()`](https://chrisycd.github.io/scDesignPop/reference/fitCopulaPop.md)
  : Fits copula for input
- [`fitMarginalPop()`](https://chrisycd.github.io/scDesignPop/reference/fitMarginalPop.md)
  : Fit marginal models for every feature
- [`fitModel()`](https://chrisycd.github.io/scDesignPop/reference/fitModel.md)
  : Fit a marginal model
- [`fitPAModel()`](https://chrisycd.github.io/scDesignPop/reference/fitPAModel.md)
  : Fit a marginal model for power analysis
- [`modifyMarginalModels()`](https://chrisycd.github.io/scDesignPop/reference/modifyMarginalModels.md)
  : Modify marginal models
- [`modifyModelPara()`](https://chrisycd.github.io/scDesignPop/reference/modifyModelPara.md)
  : Modify parameters of a glmmTMB model object
- [`normalizeSCE()`](https://chrisycd.github.io/scDesignPop/reference/normalizeSCE.md)
  : Normalize data in a SingleCellExperiment object
- [`plotCellProp()`](https://chrisycd.github.io/scDesignPop/reference/plotCellProp.md)
  : Visualize the cell type proportions across individuals
- [`plotPbulkGeno()`](https://chrisycd.github.io/scDesignPop/reference/plotPbulkGeno.md)
  : Plot pseudobulk expression vs. genotype
- [`plotReducedDimPop()`](https://chrisycd.github.io/scDesignPop/reference/plotReducedDimPop.md)
  : Dimensionality reduction and visualization for population-scale data
- [`powerAnalysis()`](https://chrisycd.github.io/scDesignPop/reference/powerAnalysis.md)
  : Perform a cell-type-specific power analysis on eQTL effects
- [`powerCICalculation()`](https://chrisycd.github.io/scDesignPop/reference/powerCICalculation.md)
  : Calculate a bootstrap confidence interval for each power
- [`runPowerAnalysis()`](https://chrisycd.github.io/scDesignPop/reference/runPowerAnalysis.md)
  : Compute power for eQTL effects in cell types
- [`simuCellProportion()`](https://chrisycd.github.io/scDesignPop/reference/simuCellProportion.md)
  : Simulate cell proportions using genotype principal components and
  population-level covariates
- [`simuNewPop()`](https://chrisycd.github.io/scDesignPop/reference/simuNewPop.md)
  : Simulate new data
- [`simulatePADesignMatrix()`](https://chrisycd.github.io/scDesignPop/reference/simulatePADesignMatrix.md)
  : Simulate new design matrix for power analysis
- [`visualizePowerCurve()`](https://chrisycd.github.io/scDesignPop/reference/visualizePowerCurve.md)
  : Visualize power as curves across study designs
- [`visualizePowerHeatmap()`](https://chrisycd.github.io/scDesignPop/reference/visualizePowerHeatmap.md)
  : Visualize power as heatmaps across study designs
