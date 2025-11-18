# m6APrediction

`m6APrediction` is an R package for predicting **N6-methyladenosine (m6A)** RNA modification sites using machine learning models such as Random Forest. 
It provides functions to encode DNA sequences, perform multiple-sample predictions, and make single-sample predictions with interpretable probability outputs.

---

## 📦 Installation

You can install the `m6APrediction` package from GitHub using the `devtools` package:

```r

devtools::install_github("zb73/m6APrediction")

```

Alternatively, you can install the package from a .tar.gz file:

```r

install.packages("path_to_your_package/m6APrediction_0.1.0.tar.gz", repos = NULL, type = "source")

```

## 🚀 Functions Overview

1. dna_encoding()
Encodes DNA sequences into a factor-based data frame representation for machine learning models.

```r

dna_encoding(c("ATCGA", "GGGTT"))

```

Description:
Each nucleotide position is encoded as a factor (A, T, C, G).
This function converts a vector of DNA strings into a data frame suitable for model training or prediction.

2. prediction_multiple()
Predicts m6A probabilities for multiple samples using a pre-trained model.

```r
ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
prediction_multiple(ml_fit, example_df, positive_threshold = 0.6)
```

## Model Performance (ROC and PRC)
## 1.ROC Curve

<img width="706" height="485" alt="ROC" src="https://github.com/user-attachments/assets/cf4d4b28-8061-4226-a311-c77d2e8005a4" />

## 2.Precision–Recall Curve

<img width="708" height="489" alt="PRC" src="https://github.com/user-attachments/assets/130a83d5-9d7f-4208-9228-3e6c2af66e75" />


## Input Columns Required:

gc_content

RNA_type (mRNA, lincRNA, lncRNA, pseudogene)

RNA_region (CDS, intron, 3'UTR, 5'UTR)

exon_length

distance_to_junction

evolutionary_conservation

DNA_5mer

## Output:
A data frame with predicted m6A probability (predicted_m6A_prob) and binary status (predicted_m6A_status).
3. prediction_single()
Predicts m6A probability for a single sample using individual input values.
```r
ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
prediction_single(
  ml_fit,
  gc_content = 0.6,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 12,
  distance_to_junction = 5,
  evolutionary_conservation = 0.8,
  DNA_5mer = "ATCGA",
  positive_threshold = 0.5
)
```
## Output:
A named vector with:

predicted_m6A_prob — probability of m6A modification

predicted_m6A_status — "Positive" or "Negative"

## 🧠 Model Requirements
This package assumes a pre-trained Random Forest model stored as an .rds file (e.g., rf_fit.rds) trained on comparable feature sets.
## 📁 Example Directory Structure
```r
m6APrediction/
├── R/
│   ├── dna_encoding.R
│   ├── prediction_multiple.R
│   ├── prediction_single.R
├── data/
│   ├── rf_fit.rds
│   └── m6A_input_example.csv
├── man/
├── DESCRIPTION
├── NAMESPACE
└── README.md
```
## 🧩 Dependencies

randomForest

base R (≥ 4.0)

## ✍️ Author
Yixuan.Pang 
(email: Yixuan.Pang23@xjtlu.edu.cn)

