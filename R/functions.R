#' @title True Positive and False Positive for SYP
#' @description This function calculates the True Positive Rate (TPR) and False Positive Rate (FPR) for the SYP dataset
#' @param data The dataset that will be used. Ensure that the columns are named appropriately (gene_name for the proteins, and SYP_ratio for the ratio of interest)
#' @param sample The name of the sample (e.g. SYP or NES) 
#' @param TP The list of true positive proteins. Ensure that the column name is gene_name.
#' @param FP The list of false positive proteins. Ensure that the column name is gene_name.
#' @param ratio The column that contains the ratio of interest (SYP or NES)
#' @param TP_flag The name of the column for the true positive flag (TP_SYP or TP_NES), so whether this protein is a TP for SYP or NES respectively
#' @return The initial dataset, including TPR and FPR 
#' @export
flag <- function(data, sample, TP, FP, gene_name, ratio, TP_flag) {
  #Flagging TP and FP proteins in the dataset
  data[[TP_flag]] <- ifelse(data[[gene_name]] %in% TP$gene_name, 1, NA)
  data$FP <- ifelse(data[[gene_name]] %in% FP$gene_name, 1, NA)
  #Calculating the median of the rations of the FP proteins
  subset <- data[!is.na(data$FP) & data$FP == 1, ]
  median <- median(subset[[ratio]], na.rm = TRUE) ##please check the name of the ratio column, it needs to match the dataset
  #Summing TP and FP
  sumTP <- sum(data[[TP_flag]], na.rm = TRUE)
  sumFP <- sum(data$FP, na.rm = TRUE)
  #Normalising and sorting
  col_name <- paste0(sample, "_FP")
  data <- data %>%
    mutate(Sample = sample,
           !!sym(col_name) := .data[[ratio]] / median) %>% ##dynamic column naming
    arrange(desc(.data[[ratio]]))
  #Initialising TPR and FPR
  data$TPR <- rep(NA, nrow(data))
  data$FPR <- rep(NA, nrow(data))
  #Calculating TPR and FPR
  x <- 1
  for (i in 1:nrow(data)) {
    if (!is.na(data[[TP_flag]][i]) && data[[TP_flag]][i] == 1) {
      data$TPR[i] <- x / sumTP
      x <- x + 1
    } else {
      data$TPR[i] <- ifelse(i > 1, data$TPR[i-1], 0)
    }
  }
  x <- 1
  for (i in 1:nrow(data)) {
    if (!is.na(data$FP[i]) && data$FP[i] == 1) {
      data$FPR[i] <- x / sumFP
      x <- x + 1
    } else {
      data$FPR[i] <- ifelse(i > 1, data$FPR[i-1], 0)
    }
  }
  return(data)
}
