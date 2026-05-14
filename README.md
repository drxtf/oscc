# OSCC Multi-omics Prognostic Signature

**Multi-omics and machine-learning framework identifies a candidate three-gene prognostic signature (PTGFRN, MICB, and IFI16) for oral squamous cell carcinoma**

[![R](https://img.shields.io/badge/R-4.3.0-blue)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.x-green)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow)](LICENSE)

## Overview

This repository contains the complete analysis code for the manuscript. The study integrates bulk transcriptomic data (TCGA-HNSC, GSE41613, GSE42743), single-cell RNA-seq (GSE172577, GSE198315), and an independent clinical cohort (n=37) to develop and validate a three-gene prognostic signature for oral squamous cell carcinoma.

## Repository Structure

```
OSCC/
├── RNA-seq/                    # Bulk transcriptomic analysis pipeline
│   ├── 01_rawdata.R            # Data acquisition (TCGA-HNSC, GEO datasets)
│   ├── 01_DEG.R                # Differential expression analysis (limma)
│   ├── 02_batch.R              # ComBat batch correction, PCA visualization
│   ├── 03_WGCNA.R              # WGCNA co-expression network analysis
│   ├── 04_Mfuzz.R              # Mfuzz temporal clustering (stages I-IV)
│   ├── 05_ML.R                 # 101-model ML benchmarking, SHAP, nomogram
│   ├── 06_oncoPredic.R         # oncoPredict drug sensitivity (GDSC2)
│   └── 07_KM.R                 # KM curves, ssGSEA, multivariate Cox, timeROC
│
├── scRNA-seq/                  # Single-cell transcriptomic analysis
│   ├── 01_rawdata.R            # scRNA-seq data loading (GSE172577, GSE198315)
│   ├── 02_celltype.R           # Integration, Harmony correction, cell annotation
│   ├── 03_copykat.R            # CopyKAT CNV inference (malignant vs normal)
│   ├── 04_Roe.R                # Ro/e enrichment, cell proportion Sankey, GSVA
│   ├── 05_marker_expression.R  # AddModuleScore & AUCell for signature genes
│   ├── 06_Pro.R                # Malignant cell subset re-analysis (Harmony)
│   ├── 07_cellchat.R           # CellChat intercellular communication (CypA pathway)
│   └── 08_pyscenic.R           # pySCENIC regulon analysis (ZEB1, IKZF1, IKZF3)
│
└── new/                        # Supplementary & validation analyses
    ├── nest.R                  # Nested cross-validation (outer 5-fold, inner 3-fold)
    ├── varPart.R               # ComBat variance decomposition analysis
    ├── valid_cross.R           # Bidirectional cross-dataset validation
    ├── tcga.R                  # TCGA-OSCC external validation (TMM, TP53 mutations)
    ├── tstage.R                # Ordinal/binary logistic regression for T-stage
    ├── tmn.R                   # TNM vs signature comparison (Cox, C-index, DCA)
    ├── significance.R          # Bootstrap CI (1000×) & PH assumption testing
    ├── relation.R              # PPIA (CypA) vs risk score Spearman correlation
    ├── calibration.R           # Calibration plots (1/3/5-year, B=1000, m=50)
    ├── immune.R                # CIBERSORT immune infiltration (22 cell types)
    └── correlation.R           # Signature genes vs CDKN2A/TP53/MKI67 correlation
```

## Pipeline Execution Order

### Phase 1: Bulk RNA-seq Discovery

```
RNA-seq/01_rawdata.R    →  Download & preprocess TCGA-HNSC, GEO datasets
RNA-seq/01_DEG.R        →  limma differential expression (|logFC|>1, FDR<0.05)
RNA-seq/02_batch.R      →  ComBat batch correction → MetaGSE.csv
RNA-seq/04_Mfuzz.R      →  Mfuzz stage-dependent clustering (10 clusters, seed=123)
RNA-seq/03_WGCNA.R      →  WGCNA (power=5, signed hybrid, minModuleSize=40)
RNA-seq/05_ML.R         →  101-model benchmarking (seed=2025, nodesize=10)
                             StepCox[forward] + RSF selected
                             SHAP interpretation, nomogram, calibration
RNA-seq/06_oncoPredic.R →  oncoPredict drug sensitivity (dasatinib, AZD7762)
RNA-seq/07_KM.R         →  KM survival, ssGSEA, multivariate Cox, timeROC
```

### Phase 2: Supplementary Analyses

```
new/varPart.R           →  ComBat efficacy: variance decomposition (paired Wilcoxon)
new/nest.R              →  Nested CV: outer K=5, inner K=3, nodesize grid=c(5,10,15)
new/valid_cross.R       →  Cross-dataset validation (GSE41613↔GSE42743)
new/immune.R            →  CIBERSORT immune deconvolution (LM22, 66 tests FDR)
new/correlation.R       →  Spearman correlation: signature vs CDKN2A/TP53/MKI67
```

### Phase 3: TCGA External Validation

```
new/tcga.R              →  TMM normalization, 3-gene surrogate score, TP53 analysis
new/tstage.R            →  Ordinal logistic regression (polr) + binary logistic (glm)
new/tmn.R               →  TNM vs signature: Cox models, timeROC, DCA, calibration
new/significance.R      →  Bootstrap AUC CI (1000 iterations) + cox.zph PH tests
new/calibration.R       →  Calibration plots (B=1000, m=50, 1/3/5-year)
new/relation.R          →  PPIA (CypA) vs risk score Spearman correlation
```

### Phase 4: Single-cell Analysis

```
scRNA-seq/01_rawdata.R          →  Load GSE172577 & GSE198315
scRNA-seq/02_celltype.R         →  QC, Harmony integration, clustering (res=0.1, dims=1:20)
scRNA-seq/03_copykat.R          →  CopyKAT CNV inference
scRNA-seq/04_Roe.R              →  Ro/e enrichment (NT→TP→TC), GSVA, dotplots
scRNA-seq/05_marker_expression.R →  Module scoring & AUCell for signature genes
scRNA-seq/06_Pro.R              →  Malignant cell re-analysis
scRNA-seq/07_cellchat.R         →  CellChat (Secreted Signaling, trim=0.1, min.cells=10)
scRNA-seq/08_pyscenic.R         →  SCENIC regulons (AUC>0.5), ZEB1/IKZF1/IKZF3 networks
```

## Software Dependencies

### R (≥4.3.0)

| Package | Version | Purpose |
|---------|---------|---------|
| limma | 3.54.0 | DEG analysis |
| sva | 3.46.0 | ComBat batch correction |
| WGCNA | 1.72.1 | Co-expression network |
| Mfuzz | 2.60.0 | Soft clustering |
| Mime | 1.0.0 | 101-model ML framework |
| randomForestSRC | ≥3.3.0 | Random Survival Forest |
| glmnet | - | LASSO regularization |
| survival | - | Cox regression, KM curves |
| survminer | - | Survival visualization |
| timeROC | - | Time-dependent ROC |
| rms | - | Calibration, nomogram |
| CIBERSORT | - | Immune deconvolution (LM22) |
| Seurat | - | scRNA-seq analysis |
| Harmony | - | Batch correction (scRNA-seq) |
| CopyKAT | - | CNV inference |
| CellChat | - | Cell-cell communication |
| oncoPredict | 0.2.1 | Drug sensitivity prediction |
| fastshap/shapviz | - | SHAP interpretation |
| ggDCA | - | Decision curve analysis |
| pwr | - | Statistical power analysis |

### Python (for pySCENIC)

| Package | Purpose |
|---------|---------|
| pySCENIC | Gene regulatory network inference |
| GRNBoost2 | GRN inference engine |
| AUCell | Regulon activity scoring |

## Key Parameters & Random Seeds

| Analysis | Parameter | Value |
|----------|-----------|-------|
| DEG | thresholds | \|logFC\| > 1, BH-adjusted p < 0.05 |
| WGCNA | power | 5 |
| WGCNA | networkType | signed hybrid |
| WGCNA | minModuleSize | 40 |
| WGCNA | mergeCutHeight | 0.15 |
| WGCNA | MAD top N genes | 10,000 |
| Mfuzz | clusters | 10 |
| Mfuzz | seed | 123 |
| ML (Mime) | train/val split | 70/30, seed=123 |
| ML (Mime) | nodesize | 10 |
| ML (Mime) | seed | 2025 |
| Nested CV | outer K | 5 |
| Nested CV | inner K | 3 |
| Nested CV | nodesize grid | {5, 10, 15} |
| Nested CV | seed | 2025 |
| Nested CV | ntree | 1000 |
| Cross-dataset | nodesize | 10 |
| Cross-dataset | seed | 2025 |
| scRNA-seq QC | nFeature_RNA | 200–7500 |
| scRNA-seq QC | MT% | < 20% |
| scRNA-seq | dimensions | 1:20 |
| scRNA-seq | resolution | 0.1 |
| CellChat | database | Secreted Signaling |
| CellChat | trim | 0.1 |
| CellChat | min.cells | 10 |
| SCENIC | AUC threshold | > 0.5 |
| Bootstrap AUC CI | iterations | 1,000 |
| Calibration | B (bootstrap) | 1,000 |
| Calibration | m (subjects per group) | 50 |
| Survival SVM | seed | 2025 |

## Key Output Files

- `MetaGSE.csv` — Batch-corrected merged GEO expression matrix (146 OSCC samples)
- `MetaGSE_OS.csv` — Overall survival annotations
- `MetaGSE_group.csv` — Dataset group labels
- `protein_cluster.txt` — Mfuzz cluster assignments
- `DEG_UP.csv` — Upregulated DEGs (227 genes)
- `101.Rdata` — Full 101-model benchmark results
- `final.ML_model.rds` — Trained StepCox[forward] + RSF model
- `TCGA_fixed_3gene_riskScore.csv` — TCGA risk scores using fixed GEO coefficients
- `CrossDataset_BestModel_Summary.csv` — Cross-dataset validation summary
- `nested_cv_full_result.rds` — Nested cross-validation results

## Data Availability

- TCGA-HNSC: https://portal.gdc.cancer.gov/projects/TCGA-HNSC
- GSE41613, GSE42743, GSE172577, GSE198315: https://www.ncbi.nlm.nih.gov/geo/
- GDSC: https://www.cancerrxgene.org/

## Citation

If you use this code, please cite the corresponding manuscript and this repository.

## License

MIT License
