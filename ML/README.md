# Machine Learning Analysis

---

## Overview

This module extends the Transcriptome-AutoPipeline framework with an exploratory machine learning analysis of retained intron (RI) splicing profiles in prostate cancer.

The analysis uses PSI (Percent Spliced In) values derived from the publicly available GEO dataset **GSE181294 (PRJNA751490)**.

Three biological groups are investigated:

- **Healthy**
- **Low-grade Tumor**
- **High-grade Tumor**

Runs annotated as single-cell in the SRA metadata were excluded to maintain consistency with the selected transcriptomic dataset.

Multiple sequencing runs belonging to the same BioSample were aggregated to avoid pseudoreplication during model evaluation.

Two supervised machine learning algorithms were implemented:

- **Logistic Regression (LR)**
- **Linear Support Vector Machine (SVM)**

Model performance was evaluated using nested stratified cross-validation, with preprocessing and hyperparameter optimization restricted to training folds.

The analysis is intended for exploratory research and does not constitute external or clinical validation.

---
## Dataset Description

The machine learning analysis was performed using the publicly available prostate cancer RNA-seq dataset **GSE181294 (PRJNA751490)**, obtained from the NCBI Gene Expression Omnibus (GEO).

The original dataset included 168 sequencing runs distributed across four biological groups: Healthy, Adjacent Normal, Low-grade Tumor, and High-grade Tumor.

For the machine learning analysis, adjacent normal samples and single-cell libraries were excluded. Only bulk RNA-seq samples from the three selected biological groups were retained.

### Final Sample Distribution

| Group | BioSamples | Sequencing Runs |
|---|---:|---:|
| Healthy | 5 | 20 |
| Low-grade Tumor | 10 | 40 |
| High-grade Tumor | 8 | 32 |
| **Total** | **23** | **92** |

Each BioSample was represented by four sequencing runs. PSI measurements were aggregated at the BioSample level to prevent data leakage caused by distributing runs from the same biological sample across training and testing folds.

The final analysis therefore included **23 independent BioSamples**, rather than treating 92 sequencing runs as independent observations.

---
## Machine Learning Workflow

The machine learning workflow was designed to investigate whether retained intron (RI) splicing profiles contain predictive information associated with prostate cancer status and tumor grade.

### 1. Data Preparation

- Import the SUPPA2-derived PSI matrix (`RI_events_with_snoRNA.tsv`).
- Integrate sequencing run metadata from `SraRunTable.csv`.
- Exclude single-cell libraries and adjacent normal samples.
- Aggregate sequencing runs belonging to the same BioSample.
- Construct a feature matrix containing RI-associated PSI values.

### 2. Feature Representation

Each biological sample is represented by a vector of PSI measurements across retained intron events.

- **Rows:** Independent BioSamples
- **Columns:** RI events
- **Features:** PSI values
- **Target:** Biological group

### 3. Machine Learning Models

Two supervised classification algorithms were investigated:

| Algorithm | Description |
|---|---|
| Logistic Regression (LR) | Linear classification baseline |
| Support Vector Machine (SVM) | Margin-based supervised classifier |
| Dummy Classifier | Reference baseline for evaluating predictive performance |

### 4. Model Validation

Model evaluation was performed using nested stratified cross-validation.

Preprocessing and model selection were performed within the training folds to reduce the risk of data leakage.

All sequencing runs from the same BioSample were aggregated before cross-validation, ensuring that the evaluation was performed at the biological sample level.

Due to the limited number of independent BioSamples, model performance should be interpreted as exploratory rather than clinically validated.

---
## Machine Learning Results

### Three-Class Classification

The primary classification task investigated whether RI-associated PSI profiles could distinguish between three biological groups:

- Healthy
- Low-grade Tumor
- High-grade Tumor

The analysis included 23 independent BioSamples after bulk RNA-seq filtering and BioSample-level aggregation.

Model performance was assessed using nested stratified cross-validation.

### Classification Performance

| Model | Accuracy | Balanced Accuracy |
|---|---:|---:|
| Logistic Regression | 0.478 | 0.483 |
| Support Vector Machine | 0.478 | 0.517 |

### Interpretation

Both models demonstrated limited predictive performance in the three-class classification task.

Although the SVM achieved slightly higher balanced accuracy than Logistic Regression, neither model provided sufficiently strong evidence for reliable discrimination between the three biological groups.

These findings suggest that the available RI-associated PSI features, combined with the limited number of independent biological samples, were insufficient to establish a robust three-class classifier.

The results should be interpreted as exploratory cross-validation estimates rather than evidence of clinical predictive utility.

---

## Limitations and Future Work

### Study Limitations

The machine learning analysis has several important limitations:

- **Limited sample size:** Only 23 independent BioSamples were available after bulk RNA-seq filtering.
- **High-dimensional feature space:** The number of RI-associated PSI features was relatively large compared with the number of biological samples.
- **Class imbalance:** The three biological groups contained unequal numbers of samples.
- **Potential batch effects:** Technical or cohort-related variation may influence PSI measurements.
- **Lack of external validation:** Model performance was estimated using cross-validation without an independent external validation cohort.

Consequently, the reported classification performance should be interpreted cautiously.

### Future Work

- [ ] **Larger cohorts** — Expand the analysis to additional independent bulk RNA-seq datasets.
- [ ] **External validation** — Evaluate trained models on independent prostate cancer cohorts.
- [ ] **Feature selection** — Investigate informative RI events using leakage-safe feature selection within cross-validation.
- [ ] **Batch effect assessment** — Examine potential technical variation across samples and datasets.
- [ ] **Biological interpretation** — Investigate relationships between predictive RI events, snoRNA loci, and associated host genes.
- [ ] **Model comparison** — Explore additional classification algorithms after establishing reliable validation procedures.

---

## Reproducibility

The machine learning module builds upon the processed outputs of the Transcriptome-AutoPipeline framework.

The analysis reuses SUPPA2-derived PSI measurements and sequencing metadata without requiring the RNA-seq preprocessing pipeline to be executed again.

### Required Input Files

| File | Description |
|---|---|
| `RI_events_with_snoRNA.tsv` | RI-associated PSI measurements generated by the transcriptomic pipeline |
| `SraRunTable.csv` | Sequencing metadata, including BioSample identifiers, library information, and biological group labels |

### Python Implementation

The machine learning implementation is available in:

`ML/gse181294_ri_ml_bulk_filtered.py`

The script provides an exploratory classification workflow based on bulk RNA-seq data, BioSample-level aggregation, and cross-validation.

The analysis distinguishes between sequencing runs and independent biological samples to reduce the risk of data leakage.

### Reproducibility Considerations

- Input datasets originate from the public GEO accession **GSE181294**.
- Single-cell libraries are excluded from the ML analysis.
- Multiple runs are aggregated at the BioSample level.
- Model evaluation uses cross-validation rather than a single random train/test split.
- The limited sample size and absence of external validation must be considered when interpreting results.

---

## Author

**Kiarash Babaei**

MSc Medical Biotechnology

University of Bologna, Italy

Research interests: Bioinformatics, Computational Genomics, Transcriptomics, Alternative Splicing, and Machine Learning in Biomedical Research.


