#' @title True Positive and False Positive Rates
#' @description This function calculates the True Positive Rate (TPR) and False Positive Rate (FPR) for the provided dataset. First, the proteomics hits are flagged as TP or FP, based on the user-provided TP and FP lists, generated based on prior annotation and literature evidence. The proteins are then normalised against the distribution of FP proteins, and ranked in descending order. Then, the TPR and FPR are calculated at each potential cutoff. The optimal cutoff will be calculated based on the ROC at a later stage.
#' @param data The dataset that will be used. 
#' @param sample The name of the sample (e.g. SYP or NES) 
#' @param TP A list or dataframe containing the true positive proteins, used for flagging. Contains proteins that should be biotinylated by TurboID. Ensure that the name of the column in the TP list matches the dataset.
#' @param FP A list or dataframe containing the false positive proteins, used for flagging. Contains proteins that do not belong in the cellular compartment being analysed. Ensure that the name of the column in the FP list matches the dataset.
#' @param gene_name The name of the column in data that contains the gene or protein names. Make sure that the name of the column in the TP and FP list matches the dataset
#' @param ratio The name of the column in data that contains the ratio of interest
#' @param TP_flag The name of the column to be created or modified in data for the true positive flag, so whether this protein is a TP for SYP or NES respectively
#' @return The initial dataset, including TPR and FPR 
#' @export
flag <- function(data, sample, TP, FP = "FP", gene_name = "gene_name", ratio, TP_flag) {
  #Flagging TP and FP proteins in the dataset
  data[[TP_flag]] <- ifelse(data[[gene_name]] %in% TP$gene_name, 1, NA)
  data$FP <- ifelse(data[[gene_name]] %in% FP$gene_name, 1, NA)
  #Calculating the median of the rations of the FP proteins
  subset <- data[!is.na(data$FP) & data$FP == 1, ]
  median <- median(subset[[ratio]], na.rm = TRUE) 
  #Summing TP and FP
  sumTP <- sum(data[[TP_flag]], na.rm = TRUE)
  sumFP <- sum(data$FP, na.rm = TRUE)
  #Normalising and sorting
  col_name <- paste0(sample, "_FP_norm")
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
#' @title ROC Curve
#' @description This function creates a ROC curve using the FPR and TPR columns in the dataset.
#' @param data The dataset that will be used. Ensure that the columns are named appropriately (FPR and TPR).
#' @param FPR The name of the column in data that contains the False Positive Rate (optional if named appropriately).
#' @param TPR The name of the column in data that contains the True Positive Rate (optional if named appropriately).
#' @param title The title of the ROC curve (optional).
#' @param col1 The color of the points in the ROC curve (optional).
#' @param col2 The color of the diagonal line in the ROC curve (optional).
#' @return A ggplot2 object showing the ROC curve.
#' @export
roc <- function(data, FPR = "FPR", TPR = "TPR", title = "ROC Curve", col1 = "#002C34", col2 = "#4F0433") {
  p <- ggplot(data = data, aes(x = FPR, y = TPR)) +
    geom_point(color = col1, shape = 19, size = 1) +
    geom_abline(intercept = 0, slope = 1, color = col2) +
    xlab("FPR") +
    ylab("TPR") +
    labs(title = title) +
    theme(plot.title = element_text(size = 20)) +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill="white"),
          axis.line = element_line("black") )
  return(p)
}
#' @title Graph to calculate the optimal Cutoff
#' @description This function creates a graph to calculate the optimal cutoff for the dataset. The optimal cutoff is the point on the ROC curve that maximises the TPR-FPR value. The function calculates the log2 of the normalised ratio and the TPR-FPR value, and plots the graph.
#' @param data The dataset that will be used. Ensure that the columns are named appropriately (FPR and TPR).
#' @param sample The name of the sample (e.g. SYP or NES)
#' @param col1 The color of the points in the graph (optional).
#' @param col2 The color of the vertical line in the graph (optional).
#' @param title The title of the graph (optional).
#' @return A ggplot2 object showing the graph to calculate the optimal cutoff.
#' @export
opt_cutoff <- function(data, sample, col1 = "#002C34", col2 = "#4F0433", title = "Optimal cutoff") {
  #Constructing column names
  norm_col <- paste0(sample, "_FP_norm")
  log_col <- paste0("log2_", sample)
  max_col <- paste0("max_", sample)
  tpr_fpr <- "TPR-FPR"
  if (!norm_col %in% names(data)) {
    stop("Column '", norm_col, "' does not exist in the dataset. Please ensure the flag() function has been run, and that the sample name matches the sample name given.")
  }
  data <- data %>%
    mutate(!!log_col := log2(as.numeric(.[[norm_col]])), 
           !!tpr_fpr := TPR - FPR)
  max_col <- max(data[[tpr_fpr]], na.rm = TRUE)
  xintersect <- data[[log_col]][data[[tpr_fpr]] == max_col]
  xintersect <- xintersect[1] #takes the first in case there are multiple. This might not be ideal and I might need to manually pick the peak instead.
  y_position <- max(data[[tpr_fpr]], na.rm = TRUE) * 0.95 #arranges where the label will be on the graph
  p <- ggplot(data, aes(x = !!sym(log_col), y = !!sym(tpr_fpr))) +
    geom_point(color = col1, shape = 19, size = 0.5) + 
    geom_vline(xintercept = xintersect, color = col2) +
    annotate("text", x = xintersect, y = y_position, label = round(xintersect, 2), 
             vjust = -0.5, color = col2, angle = 0, size = 4) +
    xlab(log_col) +
    ylab("TPR-FPR") +
    labs(x = log_col, y = tpr_fpr, title = title) +
    theme_minimal()
  
  return(p)
}
