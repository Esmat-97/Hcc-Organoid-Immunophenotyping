# HCC-Organoid-Immunophenotyping

## Transcriptomic Immunophenotyping of Hepatocellular Carcinoma Organoids

### Repository Overview

This repository contains a systems-level transcriptomic analysis pipeline for immune-associated phenotyping of hepatocellular carcinoma (HCC) patient-derived organoids (PDOs) using bulk RNA-seq data from GEO dataset GSE182593.

The project explores heterogeneous immune-related transcriptional ecosystems across HCC organoids by integrating multiple biological scoring axes including:

* Immune/Inflammatory response
* Cytotoxic activity
* Exhaustion-associated signaling
* Chemokine recruitment
* NK-associated activation
* Antigen presentation
* EMT/Stemness
* Stromal/Fibrotic remodeling

The analysis demonstrates that HCC PDOs preserve biologically distinct transcriptional states despite the absence of a complete immune microenvironment.

---

# Dataset

Dataset used:

* GEO Accession: GSE182593
* Title: Establishment of Primary Liver Cancer Organoids and Application on Sorafenib Resistance
* Organism: Homo sapiens
* Platform: Illumina HiSeq 2500

Original study focus:

* Primary liver cancer organoids
* Sorafenib resistance modeling
* Drug response characterization

This project repurposes the dataset for exploratory immunotranscriptomic systems analysis.

---

# Objectives

The main goals of this project were:

1. Infer immune-associated transcriptional states from HCC organoid RNA-seq data
2. Characterize heterogeneity across organoid samples
3. Explore relationships between:

   * Cytotoxicity
   * Exhaustion
   * Inflammatory signaling
   * EMT/stromal remodeling
4. Build an integrated immunophenotyping framework for HCC PDOs
5. Identify biologically coherent tumor ecosystem states

---

# Biological Scoring Axes

The following pathway-level scores were computed using curated marker genes:

| Score                       | Biological Interpretation                             |
| --------------------------- | ----------------------------------------------------- |
| Immune Score                | Inflamed/responder-like signaling                     |
| Exhaustion Score            | Dysfunctional/exhausted immune-associated programs    |
| Cytotoxic Score             | Killing-associated transcriptional activity           |
| Antigen Presentation Score  | Antigen processing and presentation capacity          |
| Chemokine Recruitment Score | Immune recruitment-associated signaling               |
| NK Activation Score         | Innate/NK-associated transcriptional activity         |
| Stromal Score               | Fibrotic/stromal remodeling programs                  |
| EMT Score                   | Mesenchymal transition and stemness-associated states |

---

# Key Biological Findings

## 1. Terminal Exhausted / Collapsed IFN-like State

### Cluster: orgP10

Features:

* Extremely high exhaustion
* Residual cytotoxic activity
* Minimal inflammatory signaling
* Weak chemokine recruitment
* Low NK activation

Interpretation:

* Dysfunctional cytotoxic ecosystem
* Chronic exhaustion-like phenotype
* Collapsed inflammatory recruitment program

---

## 2. EMT-Inflamed Hybrid State

### Cluster: P27RE

Features:

* Strong cytotoxic activity
* High EMT/stemness
* Moderate inflammatory signaling
* Moderate exhaustion

Interpretation:

* Hybrid mesenchymal-inflammatory phenotype
* Partial preservation of anti-tumor immune-associated programs

---

## 3. Chemokine-Dominant Inflamed State

### Cluster: P29

Features:

* High chemokine recruitment
* Strong antigen presentation
* Elevated NK signaling
* Limited cytotoxicity

Interpretation:

* Immune-recruiting but inefficient inflammatory phenotype

---

## 4. Fibrotic Mesenchymal Resistant State

### Cluster: P11

Features:

* Extremely high stromal/fibrotic signaling
* High EMT
* High exhaustion
* Weak cytotoxicity

Interpretation:

* Fibrotic immunosuppressive ecosystem
* Mesenchymal resistance-associated phenotype

---

## 5. Mixed Inflamed-Stromal State

### Cluster: p26RE

Features:

* Highest immune/responder-like scores
* High stromal remodeling
* High EMT
* Moderate cytotoxicity

Interpretation:

* Inflamed but immunorestrained ecosystem

---

## 6. Innate-like Immune-Cold State

### Clusters: P12 / P16RE

Features:

* Strong NK-associated activity
* Weak adaptive inflammatory signaling
* Low cytotoxicity

Interpretation:

* Innate-like activation without coordinated adaptive immune programs

---

# Major Systems-Level Conclusions

## Immune activation is multidimensional

The analysis revealed partial uncoupling between:

* Cytotoxicity
* Inflammatory signaling
* NK activation
* EMT
* Stromal remodeling

This demonstrates that HCC organoids preserve multiple distinct transcriptional immune-associated ecosystems rather than a single linear immune continuum.

---

# Visualizations Included

The repository includes:

* Barplots for all biological scores
* Integrated clustered heatmaps
* Immune state comparison plots
* Multi-axis immunophenotyping visualizations
* Hierarchical clustering analyses

---

# Technologies Used

* Python
* Pandas
* NumPy
* Seaborn
* Matplotlib
* SciPy


---

# Scientific Scope

This project should be interpreted as:

## Exploratory Transcriptomic Immunophenotyping

The inferred immune-associated states are:

* transcriptomic
* computational
* hypothesis-generating

and do not represent:

* experimentally validated immune infiltration
* confirmed checkpoint inhibitor response
* direct immune-cell quantification

---

# Potential Future Directions

* Single-cell validation
* Immune deconvolution approaches
* GSVA/ssGSEA pathway enrichment
* Survival association studies
* Drug-response integration
* Multi-omics expansion
* Spatial transcriptomics validation

---

# Citation

If you use this repository or analysis framework, please cite the original GEO dataset and associated publication:

Chen X, Wu T, Xian L, Ma L et al.
"circGLS2 inhibits hepatocellular carcinoma recurrence via regulating hsa-miR-222-3p-PTEN-AKT signaling."
Signal Transduction and Targeted Therapy (2023)

---

# Author

Mohamed Esmat

Graduate Researcher in Cancer Biology & Immunology

---

# License

MIT License
