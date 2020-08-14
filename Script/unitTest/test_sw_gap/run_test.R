library(testthat)

suppressPackageStartupMessages(c(library(Biostrings),
                               library(BiocGenerics),
                               library(parallel),
                               library(S4Vectors),
                               library(seqinr),
                               library(stringr),
                               library(stringi)))

# file = "~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script/sw_gap.R"
main = function(file,  dir){
  
  Window <<- 6
  Wall   <<- 12
  
  source(file)
  test_results = test_dir(dir, reporter = "summary")
  test_results
}

args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2])



