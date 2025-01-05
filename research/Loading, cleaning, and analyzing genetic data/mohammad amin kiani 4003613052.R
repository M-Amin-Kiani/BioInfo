library(ggplot2)
library(reshape2)


# Set the working directory
setwd("C:/Users/Lenovo/Desktop/Train1")
print(getwd())

# Try to read the CSV file
data_file_path <- "data.csv"
geneExp <- read.csv(data_file_path, header = FALSE)

# Check if the file was successfully read
if (length(geneExp) > 0) {
  print(is.data.frame(geneExp))
  print(dim(geneExp))
  
  # Check for non-numeric cells
  non_numeric <- apply(geneExp, 2, function(x) any(is.na(as.numeric(x))))
  
  # Display columns that contain non-numeric values
  if (any(non_numeric)) {
    cat("Non-numeric columns found:", paste(names(geneExp)[non_numeric], collapse = ", "), "\n")
  } else {
    cat("All columns are numeric.\n")
  }
  
  # Replace "#REF!" values with 0
  geneExp[geneExp == "#REF!"] <- 0
  
  # Extract just the numeric data into a matrix
  rownames(geneExp) <- geneExp$V1
  geneExp_matrix <- as.matrix(geneExp[, 2:16])
  
  # Apply numeric conversion and whitespace trimming
  numeric_matrix <- apply(geneExp_matrix, c(1, 2), function(x) as.numeric(trimws(x)))
  
  # Create heatmaps
  heatmap(numeric_matrix)
  log2_matrix <- log2(numeric_matrix + 1)
  heatmap(log2_matrix)
} else {
  print("Error: Unable to read the CSV file.")
}
