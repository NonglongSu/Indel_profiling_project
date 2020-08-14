#!/bin/bash

#Variables
RSCRIPT=Rscript --vanilla 
RM=rm -i -f

Dir=Raw_data.2.outgroups


#default:all
#all:    mafft_upc prank_upc mafft_map prank_map mafft_plus prank_plus mafft_phase prank_phase mafft_dis prank_dis mafft_edge prank_edge 

.PHONY:  all default clean 
.PHONY:  mafft_upc prank_upc mafft_map prank_map mafft_plus prank_plus mafft_phase prank_phase mafft_dis prank_dis mafft_edge prank_edge  


#Check ../../direct.0.sh to build up required directories. 

#######################################################

NAMELIST:=$(shell cat ../../$(Dir)/nameList.1.txt )


Mafft/updated_cds:	 mafft_upc
Prank/updated_cds:	 prank_upc
Mafft/mapped_cds_Out:	 mafft_map
Prank/mapped_cds_Out:	 prank_map
Mafft/mapped_cds_plus:   mafft_plus
Prank/mapped_cds_plus: 	 prank_plus


mafft_upc:	 $(patsubst %, Mafft/updated_cds/%.fa,    $(basename $(NAMELIST)))
prank_upc:	 $(patsubst %, Prank/updated_cds/%.fa,    $(basename $(NAMELIST)))
mafft_map:	 $(patsubst %, Mafft/mapped_cds_Out/%.fa, $(basename $(NAMELIST)))
prank_map:	 $(patsubst %, Prank/mapped_cds_Out/%.fa, $(basename $(NAMELIST)))

##########################################################


#Update the gaps
Mafft/updated_cds/%.fa: ../../Script/update_gap.R   ../../$(Dir)/mapped_cds_mafft/%.fa Mafft/updated_cds/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)
.SECONDARY:Mafft/updated_cds/%.fa

Prank/updated_cds/%.fa: ../../Script/update_gap.R   ../../$(Dir)/mapped_cds_prank/%.fa Prank/updated_cds/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)
.SECONDARY:Prank/updated_cds/%.fa


#Find the alternative alignments
Mafft/mapped_cds_Out/%.fa: ../../Script/pseudo_seq.R   Mafft/updated_cds/%.fa Mafft/mapped_cds_Out/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)
.SECONDARY:Mafft/mapped_cds_Out/%.fa

Prank/mapped_cds_Out/%.fa: ../../Script/pseudo_seq.R   Prank/updated_cds/%.fa Prank/mapped_cds_Out/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)
.SECONDARY:Prank/mapped_cds_Out/%.fa


#############################################
#Generate nameFilter.txt and mapped_cds_In_Out/
#bash test_human_mouse_rat/mafft.nameFilter.sh 
#bash test_human_mouse_rat/prank.nameFilter.sh 	
NAMELIST.1:=$(shell cat Mafft/nameFilter.txt) 
NAMELIST.2:=$(shell cat Prank/nameFilter.txt)

mafft_plus:	 $(patsubst %, Mafft/mapped_cds_plus/%.fa,$(basename $(NAMELIST.1)))
prank_plus:	 $(patsubst %, Prank/mapped_cds_plus/%.fa,$(basename $(NAMELIST.2)))
	
############################################### 

#Add +++ to mapped_cds_In_Out/ and redirect to mapped_cds_plus/
Mafft/mapped_cds_plus/%.fa: ../../Script/update_gap_plus.R  Mafft/mapped_cds_Out/%.fa Mafft/mapped_cds_In_Out/%.fa Mafft/mapped_cds_plus/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^)
.SECONDARY:Mafft/mapped_cds_plus/%.fa

Prank/mapped_cds_plus/%.fa: ../../Script/update_gap_plus.R  Prank/mapped_cds_Out/%.fa Prank/mapped_cds_In_Out/%.fa Prank/mapped_cds_plus/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^)
.SECONDARY:Prank/mapped_cds_plus/%.fa

#########################################################################
	 



clean:
	
	@$(RM) Mafft/updated_cds/*
	@$(RM) Mafft/mapped_cds_Out/*
	@$(RM) Mafft/mapped_cds_plus/*
	@$(RM) Mafft/mapped_cds_In_Out/*
	@$(RM) Prank/updated_cds/*
	@$(RM) Prank/mapped_cds_out/*
	@$(RM) Prank/updated_cds_plus/*
	@$(RM) Prank/updated_cds_In_Out/*

