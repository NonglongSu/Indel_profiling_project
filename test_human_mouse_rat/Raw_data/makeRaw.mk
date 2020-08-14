#!/bin/bash

#Set var
RSCRIPT = Rscript --vanilla
MAFFT   = mafft --maxiterate 1000 --globalpair
PRANK   = prank -F -f=fasta
RM    = rm -i -f
S     = ../../Script

DATABASE   := geneId.txt
NAMELIST   := $(shell test -f $(DATABASE) && cat  $(DATABASE) | cut -f 1 | sed -e 's/"//g')
NAMELIST1  := $(shell test -f nameList.txt && cat nameList.txt)

#########################
default: all
all:     

.PHONY: default all clean \
        homoCall cds test_mul3 filter test_N test_stop test_preStop \
	addStop aa aaFix guideTree mafft prank mapped_mafft mapped_prank

###########################################PART I

homoCall: $(S)/$(HOMO) $(INPUT)
	$(RSCRIPT) $< $(word 2,$^) geneId.txt   

cds_seq: cds
cds: $(patsubst %, cds_seq/%.fa, $(NAMELIST))
cds_seq/%.fa: $(S)/cds_get.R $(DATABASE) 
	$(RSCRIPT) $< $@ $(word 2,$^) cds_seq/ 
.SECONDARY: cds_seq/%.fa

###########################################Quality control

#Record fasta files which length are not-multiple-of-3.
test_mul3: $(S)/test_mul_3.R cds_seq/
	$(RSCRIPT) $< $(word 2,$^) QC/multi_no_3.txt

#Remove them and update the database.
filter: $(S)/filter.sh cds_seq/ 
	bash $< $(word 2,$^)

#Record files that have 'N' fuzzy code.
test_N: $(S)/test_NNN.R cds_seq/
	$(RSCRIPT) $< $(word 2,$^) QC/NNN.txt

# Record files without stop codons.
test_stop: $(S)/test_std.R cds_seq/
	$(RSCRIPT) $< $(word 2,$^) QC/noStop.txt

#Record files with pre stop codons.
test_preStop: $(S)/test_preStop.R cds_seq/
	$(RSCRIPT) $< $(word 2,$^) QC/preStop.txt 


########################################################PART II

aa_seq:				aa	
Alignments/mafft_out:		mafft
Alignments/prank_out:		prank
mapped_cds_mafft:		mapped_mafft
mapped_cds_prank:		mapped_prank


aa:      	$(patsubst %, aa_seq/%.fa, $(NAMELIST1))
mafft:   	$(patsubst %, Alignments/mafft_out/%.fa, $(NAMELIST1))
prank:   	$(patsubst %, Alignments/prank_out/%.fa, $(NAMELIST1))
mapped_mafft:   $(patsubst %, mapped_cds_mafft/%.fa, $(NAMELIST1))
mapped_prank:   $(patsubst %, mapped_cds_prank/%.fa, $(NAMELIST1))

aa_seq/%.fa:         		 | cds_seq/%.fa
Alignments/mafft_out/%.fa:       | aa_seq/%.fa      
Alignments/prank_out/%.fa:	 | aa_seq/%.fa

mapped_cds_mafft/%.fa:	 | Alignments/mafft_out/%.fa cds_seq/%.fa
mapped_cds_prank/%.fa:	 | Alignments/prank_out/%.fa cds_seq/%.fa


#Add stop codon.
addStop: $(S)/add_stop.R QC/noStop.txt cds_seq/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^)

#Translation           
aa_seq/%.fa:$(S)/dna2aa.R cds_seq/%.fa 
	$(RSCRIPT) $< $(word 2,$^) aa_seq/ ${ARRAY}          
.SECONDARY: aa_seq/%.fa

#Convert '*' to 'X' in pre stop codon.
aaFix: $(S)/transform_preStop.R QC/preStop.txt aa_seq/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) 


#MAFFT
Alignments/mafft_out/%.fa: aa_seq/%.fa
	$(MAFFT) $< > $@
.SECONDARY: Alignments/mafft_out/%.fa

#PRANK
guideTree: $(S)/make_guide_tree.R $(INPUT)
	$(RSCRIPT) $< $(word 2,$^) > guide.tree

Alignments/prank_out/%.fa: aa_seq/%.fa                     
	$(PRANK) -d=$< -t=guide.tree -o=Alignments/prank_out/$*
	test -f Alignments/prank_out/$*.best.fas && \
	samtools faidx Alignments/prank_out/$*.best.fas $(shell cut -f3 $(INPUT)) > $@ && \
	rm Alignments/prank_out/$*.best*

.SECONDARY: Alignments/prank_out/%.fa


#Abnormal mapping
#Do not forget to change the file.path in mapping_cds2aa.R

mapped_cds_mafft/%.fa: $(S)/map_cds2aa.R  Alignments/mafft_out/%.fa cds_seq/%.fa 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) mapped_cds_mafft/
.SECONDARY: mapped_cds_mafft/%.fa

mapped_cds_prank/%.fa: $(S)/map_cds2aa.R  Alignments/prank_out/%.fa cds_seq/%.fa 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) mapped_cds_prank/
.SECONDARY: mapped_cds_prank/%.fa


#Run [bash nameCheck.sh emptyCheck.sh] to fix any broken pipe.
##############################################################PART III


clean:
	@$(RM) cds_seq/*
	@$(RM) aa_seq/*
	@$(RM) Alignments/mafft_out/*
	@$(RM) Alignments/prank_out/*
	@$(RM) mapped_cds_mafft/*
	@$(RM) mapped_cds_prank/*

help:
	@echo "help me, master Ziqi!"	






