#Test each function in pseudo_seq.R
context("Pseudo_Seq.R script testing process start:")

##Pseduo data 
seq1 = "AAT===AAACAAAGAATGCTTACTGT---ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
seq2 = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTG---ATCGCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"

#global vars
wall <<- 6
lwall <<- wall -1
rwall <<- wall +3

##Predicted outcomes
#index
res1 = c(26,122)
res2 = 107

#alignemnt score (Do not use circulating decimal)
res3 = 16.6666667
res4 = 41.6666667

#left/right slide mode
seq1.l = "AAT===AAACAAAGAATGCTTACTG---TATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
seq2.l = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTT---GATCGCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"
  
seq1.r = "AAT===AAACAAAGAATGCTTACTGTA---TAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
seq2.r = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTGA---TCGCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"

#Best alignment mode
#unchanged
best_aligned.m = "AAT===AAACAAAGAATGCTTACTGT+++ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCA---TAATAGGGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
#changed
best_aligned.m2= "AAT===AAACAAAGAATGCTTACTGT---ATAAGGCTTACTGTTCTAGCG===ATCACCGCG===TCATGTCTAGTTATGAACGGC------GGTTTAACATTGAATAGCAAGGCACTTCCATAATAG---GGCCGTC===GTAATTGTCTAATATAG------ATAGTA==="
best_aligned.r = "TAA------AA===AATTTGATGCTACATTGGATGAGTCTACTTCGAGCGCGCCGCATCGATTGCAAGAGCAGTGTTGCCT===AAGAGCCGTTAGATGCGTCGTTGATC---GCGTCCGATAATTCGGGAGTTGTGC===CCCAATATTTAATATGATGA===TAGCTATAA"

# Test Index func.
test_that("Test Index() ",{
  expect_equal(Index(seq1),res1)
  expect_equal(Index(seq2),res2)
})

# Test Align func. 
wid.1     = substr(seq1, start =res1[1]-lwall, stop = res1[1]+rwall)
wid_ref.1 = substr(seq2, start =res1[1]-lwall, stop = res1[1]+rwall)
wid.2     = substr(seq2, start =res2[1]-lwall, stop = res2[1]+rwall)
wid_ref.2 = substr(seq1, start =res2[1]-lwall, stop = res2[1]+rwall)

test_that("Test Align() ",{
  expect_equal(Align(wid.1,wid_ref.1),res3)
  expect_equal(Align(wid.2,wid_ref.2),res4)
})

# Test left_slide func.
test_that("Test left_slide() ",{
  expect_equal(left_slide(seq1,str_convert(seq1),res1[1]),seq1.l)
  expect_equal(left_slide(seq2,str_convert(seq2),res2[1]),seq2.l)
})

# Test right_slide func.
test_that("Test right_slide() ",{
  expect_equal(right_slide(seq1,str_convert(seq1),res1[1]),seq1.r)
  expect_equal(right_slide(seq2,str_convert(seq2),res2[1]),seq2.r)
})

# Test Merge func.
test_that("Test Merge()",{
   expect_equal(Merge(res1[1],seq1,seq2),best_aligned.m)
   expect_equal(Merge(res1[2],seq1,seq2),best_aligned.m2)
   expect_equal(Merge(res2[1],seq2,seq1),best_aligned.r)
})





