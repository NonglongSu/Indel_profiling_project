#!/bin/bash

default:all
all:    mafft_upc mafft_map mafft_plus

.PHONY: all default clean 
.PHONY: mafft_upc mafft_map mafft_plus 


#Set up vars
RSCRIPT = Rscript --vanilla 
RM = rm -i -f

scr = ../../Script
Dir = ../Raw_data

#####################################################PART I

NAMELIST := $(shell cat $(Dir)/nameList.txt)


#Remember to switch the bash file between Mafft and Mafft_Mul        >>>>>>>>>>>>>>
#NAMELIST.1 := $(shell bash ../mafft.nameFilter.sh) 
#NAMELIST.1 := $(shell cat Mafft/nameFilter.txt)

#NAMELIST.2 := $(shell bash ../prank.nameFilter.sh)
#NAMELIST.2 := $(shell cat Prank/nameFilter.txt)


mafft_upc:	 $(patsubst %, Mafft/updated_cds/%.fa,  $(basename $(NAMELIST)))
mafft_map:	 $(patsubst %, Mafft/mapped_cds/%.fa,   $(basename $(NAMELIST)))
#mafft_plus:	 $(patsubst %, Mafft/mapped_cds_plus/%.fa,$(basename $(NAMELIST.1)))

Mafft/updated_cds/%.fa:  | $(Dir)/mapped_cds_mafft/%.fa
Mafft/mapped_cds/%.fa:   | Mafft/updated_cds/%.fa 


#prank_upc:	 $(patsubst %, Prank/updated_cds/%.fa,    $(basename $(NAMELIST)))
#prank_map:	 $(patsubst %, Prank/mapped_cds/%.fa,     $(basename $(NAMELIST.2)))
#prank_plus:	 $(patsubst %, Prank/mapped_cds_plus/%.fa,$(basename $(NAMELIST.2)))

#constant denominator                                                >>>>>>>>>>>>>>>>>>      
#mafft_map_cD:	$(patsubst %, Mafft/mapped_cds_out_cD/%.fa, $(basename $(NAMELIST)))
#mafft_plus_cD: $(patsubst %, Mafft/mapped_cds_plus_cD/%.fa,$(basename $(NAMELIST.1)))

#######################################################PART II

#Update all legal gaps
Mafft/updated_cds/%.fa: $(scr)/update_gap.R $(Dir)/mapped_cds_mafft/%.fa Mafft/updated_cds/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:Mafft/updated_cds/%.fa

#Find the best phase
Mafft/mapped_cds/%.fa:  $(scr)/sw_gap.R Mafft/updated_cds/%.fa Mafft/mapped_cds/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:Mafft/mapped_cds/%.fa


#Prank/updated_cds/%.fa: $(scr)/update_gap.R $(Dir)/mapped_cds_prank/%.fa Prank/updated_cds/ 
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
#.SECONDARY:Prank/updated_cds/%.fa
#
#Prank/mapped_cds/%.fa:  $(scr)/pseudo_seq.R Prank/updated_cds/%.fa Prank/mapped_cds_Out/
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
#.SECONDARY:Prank/mapped_cds/%.fa


#############################################################
#Constant denominator
#
#Mafft/mapped_cds_Out_cD/%.fa: $(scr)/Script/pseudo_seq_cD.R Mafft/updated_cds_cD/%.fa Mafft/mapped_cds_Out_cD/ 
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)
#.SECONDARY:Mafft/mapped_cds_Out_cD/%.fa
#
#############################################################

#Add +++ to mapped_cds_In_Out/ and redirect to mapped_cds_plus/
Mafft/mapped_cds_plus/%.fa: ../../Script/update_gap_plus.R  Mafft/mapped_cds_Out/%.fa Mafft/mapped_cds_In_Out/%.fa Mafft/mapped_cds_plus/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^)
.SECONDARY:Mafft/mapped_cds_plus/%.fa

Prank/mapped_cds_plus/%.fa: ../../Script/update_gap_plus.R  Prank/mapped_cds_Out/%.fa Prank/mapped_cds_In_Out/%.fa Prank/mapped_cds_plus/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^)
.SECONDARY:Prank/mapped_cds_plus/%.fa

#########################################################################
#Create this folder for calulating the displacement distribution later. 
#Add +++ to mapped_cds_In_Out_cD/ and redirect to mapped_cds_plus_cD/
Mafft/mapped_cds_plus_cD/%.fa: ../../Script/update_gap_plus.R  Mafft/mapped_cds_Out_cD/%.fa Mafft/mapped_cds_In_Out_cD/%.fa Mafft/mapped_cds_plus_cD/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^)
.SECONDARY:Mafft/mapped_cds_plus_cD/%.fa




clean:
	
	@$(RM) Mafft/updated_cds/*
	@$(RM) Mafft/mapped_cds/*
	@$(RM) Prank/updated_cds/*
	@$(RM) Prank/mapped_cds/*

	
