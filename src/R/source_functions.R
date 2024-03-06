# Source functions called by process_input_data.R

#' Write CLS file
#' 
#' Function to write a GSEA phenotype labels .cls file.
#' @param groups Dataframe defining sample_id and group label.
#' @param out_basename Path and base filename for output .cls file
#' @returns No return value, writes to file and prints a success message.
write_cls <- function(groups, out_basename) {
  line1_a<- nrow(groups)
  line1_b <- length(unique(groups$group))
  line1 <- paste(c(line1_a, line1_b, 1), collapse = " ")
  line2_a <- paste(levels(groups$group), collapse = " ")
  line2 <- paste(c("#", line2_a), collapse = " ")
  line3 <- paste(groups$group, collapse = " ")
  cls_file <- str_glue("{out_basename}_phenolabels.cls") 
  
  write(x = line1, file = cls_file, ncolumns = 3)
  write(x = line2, file = cls_file, ncolumns = 3, append = TRUE)
  write(x = line3, file = cls_file, ncolumns = nrow(groups), append = TRUE)
  
  print(paste("Saved phenotype labels to:", cls_file))
}

# # function documentation (run once)
# setwd("~/Documents/jobs/code_sample_projects/GSEAforR/src")
# cat("Package: GSEAforR\n", file = "DESCRIPTION")
# cat("Version: 1.0.0\n", file = "DESCRIPTION", append = TRUE)
# roxygen2::roxygenise()
