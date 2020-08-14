#!/bin/bash

#Variables
RSCRIPT=Rscript --vanilla 
RM=rm -i -f
window=6
Dir=Raw_data.2.outgroups


#default:all
#all:    mafft_upc prank_upc mafft_map prank_map mafft_anc prank_anc  

.PHONY:  all default clean 
.PHONY:  mafft_upc prank_upc mafft_map prank_map mafft_anc prank_anc 


#######################################################

NAMELIST:=$(shell cat ../../$(Dir)/nameList.1.txt )

#NAMELIST.1:=$(shell bash ../../mafft.nameFilter.sh) 
#NAMELIST.2:=$(shell bash ../../prank.nameFilter.sh)


mafft_upc:	 $(patsubst %, Mafft_Mul/updated_cds/%.fa,    $(basename $(NAMELIST)))
prank_upc:	 $(patsubst %, Prank_Mul/updated_cds/%.fa,    $(basename $(NAMELIST)))
mafft_map:	 $(patsubst %, Mafft_Mul/mapped_cds_Out/%.fa, $(basename $(NAMELIST)))
prank_map:	 $(patsubst %, Prank_Mul/mapped_cds_Out/%.fa, $(basename $(NAMELIST)))

#Generate nameFilter.txt 
#bash ../../mafft.nameFilter.sh
#bash ../../prank.nameFilter.sh
NAMELIST.1:=$(shell cat Mafft_Mul/nameFilter.txt)
#NAMELIST.2:=$(shell cat Prank_Mul/nameFilter.txt)

mafft_anc:	 $(patsubst %, Mafft_Mul/mapped_anc/%.fa,$(basename $(NAMELIST.1)))
prank_anc:	 $(patsubst %, Prank_Mul/mapped_anc/%.fa,$(basename $(NAMELIST.2)))

##########################################################


#Update the gaps
Mafft_Mul/updated_cds/%.fa: ../../Script/update_gap_mul.R   ../../$(Dir)/mapped_cds_mafft_Tmp/%.fa Mafft_Mul/updated_cds/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)
.SECONDARY:Mafft_Mul/updated_cds/%.fa

Prank_Mul/updated_cds/%.fa: ../../Script/update_gap_mul.R   ../../$(Dir)/mapped_cds_prank_Tmp/%.fa Prank_Mul/updated_cds/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)
.SECONDARY:Prank_Mul/updated_cds/%.fa


#Find the best alignment
Mafft_Mul/mapped_cds_Out/%.fa: ../../Script/pseudo_seq_mul.R   Mafft_Mul/updated_cds/%.fa Mafft_Mul/mapped_cds_Out/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)
.SECONDARY:Mafft_Mul/mapped_cds_Out/%.fa

Prank_Mul/mapped_cds_Out/%.fa: ../../Script/pseudo_seq_mul.R   Prank_Mul/updated_cds/%.fa Prank_Mul/mapped_cds_Out/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)
.SECONDARY:Prank_Mul/mapped_cds_Out/%.fa


#Find the common ancestor and update the focal sequences
Mafft_Mul/mapped_anc/%.fa: ../../Script/add_com_anc.R  Mafft_Mul/mapped_cds_Out/%.fa Mafft_Mul/mapped_anc/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) 
.SECONDARY:Mafft_Mul/mapped_anc/%.fa

Prank_Mul/mapped_anc/%.fa: ../../Script/add_com_anc.R  Prank_Mul/mapped_cds_Out/%.fa Prank_Mul/mapped_anc/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^)
.SECONDARY:Prank_Mul/mapped_anc/%.fa

#########################################################################




clean:
	
	@$(RM) Mafft_Mul/updated_cds/*
	@$(RM) Mafft_Mul/mapped_cds_Out/*
	@$(RM) Prank_Mul/updated_cds/*
	@$(RM) Prank_Mul/mapped_cds_out/*

