# Profliling of Indel phases in coding regions
* Find the **real indels ** that cannot be determined due to the limitation of current software.  
* Use indel phases  

## Requirements
* R 3.4.4 + 
* Bioconductor 3.6 +
* Ensembl Genes 100 (Human Genes GRCh38.p13)
* Mafft --version V.7.407
  Prank --version v.170427  

## Workflow
## Part I
### 1. Find homologs in multiple species (Ex. human/mouse/rat) 
make homoCall 
* human/mouse/rat &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;-- 14547
* human/hamster/mouse/rat -- 13890
* mouse/human/macaque &nbsp; &nbsp; &nbsp;--14422

### 2. Extract all protein coding sequences in each species.
make cds  
make filter  
make addStop  
* human/mouse/rat &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;-- 14532  
* human/hamster/mouse/rat -- 13881
* mouse/human/macaque &nbsp; &nbsp; &nbsp;--14418

### 3. Quality control test
make test_mul_3  
make test_N  
make test_stop  
make test_preStop  

### 4. Translation of nucleotide to amino acid for each sequence.
make aa  
make aaFix  

### 5. Multiple sequence alignment of cds.
make mafft  
make prank  

### 6. Abnormal mapping (From aa --> dna)
make mapped_mafft  
make mapped_prank  

## Part II
### 7. Include gaps of length = {3, 6, 9, 12} only.
make mafft_upc   
make prank_upc  

### 8. Find the "true" gaps based on sliding-window method. 
make mafft_map  
* h/m/r >>>>
* Data_3   : 2034
* Data_6   : 1712 
* Data_9   : 1652
* Data_12  : 1605

* h/c/m/r >>>>

* m/h/m   >>>>

make prank_map




## Figure
### 1. Count the proportion of phase 0, 1, 2 indels of focal species together.
make mafft_phase
make prank_phase

### 2. Generate displacement bar plot between prior-post indels.
make mafft_dis
make prank_dis

### 3. Generate edge-bias HTML report
make mafft_edge_report
make prank_edge_report

### 4. Count the proportion of phase 0, 1, 2 indels of focal species seperately.
make mafft_phase_plot
make prank_phase_plot

### 5. Count the effective proportion of phase 1 & 2 indels.
make_eff_phase
