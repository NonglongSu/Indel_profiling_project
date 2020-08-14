library(testthat)

suppressPackageStartupMessages(c(library(Biostrings),
                                 library(BiocGenerics),
                                 library(parallel),
                                 library(S4Vectors),
                                 library(seqinr),
                                 library(stringi),
                                 library(XVector)))

# file = "~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script/update_gap.R"
main = function(file, dir){
 
  Window <<- 6
  Wall   <<- 12
  
  source(file)
  test_results = test_dir(dir, reporter = "summary")
  
  #Print in terminal
  test_results
}

args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])

