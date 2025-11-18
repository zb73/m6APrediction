#' DNA Sequence Encoding
#'
#' This function encodes DNA strings as a data frame where each nucleotide is represented as a factor.
#' Each position in the sequence is labeled as a factor with levels corresponding to the nucleotides A, T, C, and G.
#'
#' @param dna_strings A character vector of DNA sequences to be encoded.
#' @return A data frame where each column represents a position in the DNA sequence and each element is a factor
#'         with levels "A", "T", "C", "G" corresponding to the nucleotides.
#' @examples
#' \dontrun{
#' dna_encoding(c("ATCGA", "GGGTT"))
#' }
#' @import randomForest
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Multiple Sample Prediction for m6A Sites
#'
#' This function takes a trained machine learning model and a data frame of features and predicts m6A probabilities
#' for multiple samples. It uses a pre-trained model to infer probabilities for each sample, assigns a status of
#' "Positive" or "Negative" based on a given threshold, and returns the input data frame with predicted probabilities
#' and status.
#'
#' @param ml_fit A trained machine learning model (e.g., Random Forest) for predicting m6A sites.
#' @param feature_df A data frame containing the features for prediction. Must include columns:
#'        "gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
#'        "evolutionary_conservation", and "DNA_5mer".
#' @param positive_threshold A numeric value between 0 and 1, used to classify the predicted m6A status as "Positive"
#'        if the probability is greater than the threshold (default is 0.5).
#' @return A data frame with the same columns as `feature_df` plus two new columns:
#'         "predicted_m6A_prob" (the predicted probability of m6A) and "predicted_m6A_status"
#'         (either "Positive" or "Negative" based on the threshold).
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package="m6APrediction"))
#' prediction_multiple(ml_fit, example_df, positive_threshold = 0.6)
#' @export
#' @import randomForest
#' @importFrom stats predict
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Single Sample Prediction for m6A Sites
#'
#' This function takes individual feature values and uses a pre-trained machine learning model to predict the
#' m6A probability and status for a single sample. It calls the `prediction_multiple` function with a one-row data frame
#' and returns the predicted values as a named vector.
#'
#' @param ml_fit A trained machine learning model (e.g., Random Forest) for predicting m6A sites.
#' @param gc_content A numeric value for the GC content of the sample.
#' @param RNA_type A character string specifying the RNA type. One of: "mRNA", "lincRNA", "lncRNA", "pseudogene".
#' @param RNA_region A character string specifying the RNA region. One of: "CDS", "intron", "3'UTR", "5'UTR".
#' @param exon_length A numeric value for the length of the exon in the sample.
#' @param distance_to_junction A numeric value for the distance to the nearest junction.
#' @param evolutionary_conservation A numeric value for the evolutionary conservation score of the region.
#' @param DNA_5mer A character string representing the 5-mer DNA sequence.
#' @param positive_threshold A numeric value between 0 and 1, used to classify the predicted m6A status as "Positive"
#'        if the probability is greater than the threshold (default is 0.5).
#' @return A named vector with two elements: "predicted_m6A_prob" (the predicted probability of m6A)
#'         and "predicted_m6A_status" (either "Positive" or "Negative" based on the threshold).
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
#'                   exon_length = 12, distance_to_junction = 5, evolutionary_conservation = 0.8,
#'                   DNA_5mer = "ATCGAT", positive_threshold = 0.5)
#' @export
#' @import randomForest
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}

