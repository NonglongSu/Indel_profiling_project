library(testthat)

suppressPackageStartupMessages(c(library(Biostrings),
                                 library(BiocGenerics),
                                 library(parallel),
                                 library(S4Vectors),
                                 library(seqinr),
                                 library(stringr)))

file = "~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script/update_gap.R"

source(file)
test_results = test_dir("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script/Test_update", reporter="summary")
test_results
