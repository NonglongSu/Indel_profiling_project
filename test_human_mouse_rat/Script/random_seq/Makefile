#!/bin/bash

#Create a simulated dna sequence for alignment

default:all

all: md extra_nog sim_align rm_gap aa aa_fix mafft mapped_mafft mafft_upc mafft_map new mafft_plus

.PHONY: md extra_nog sim_align rm_gap aa aa_fix mafft mapped_mafft mafft_upc mafft_map new mafft_plus mafft_dis 
.PHONY: clean

################################################

RSCRIPT=Rscript --vanilla 
MAFFT=mafft --maxiterate 1000 --preservecase --globalpair
array="c('mouse','rat')"

Dir=../../Raw_data.2.outgroups/mafft_sim
NAMELIST:=$(shell head -2000 nogaps.csv)
window = 6

sim_align:	$(patsubst %, $(Dir)/mafft_ali/%.fa,      $(basename $(NAMELIST)))
aa:             $(patsubst %, $(Dir)/aa_seq/%.fa,         $(basename $(NAMELIST)))
aa_fix:         $(patsubst %, $(Dir)/aa_seq_fix/%.fa,     $(basename $(NAMELIST)))
mafft:          $(patsubst %, $(Dir)/mafft_reali/%.fa,    $(basename $(NAMELIST)))
mapped_mafft:   $(patsubst %, $(Dir)/mafft_mapped/%.fa,   $(basename $(NAMELIST)))
mafft_upc:      $(patsubst %, $(Dir)/updated_gaps/%.fa,   $(basename $(NAMELIST)))
mafft_map:      $(patsubst %, $(Dir)/sw_gaps/%.fa, 	  $(basename $(NAMELIST)))

NAMELIST.1:=$(shell cat $(Dir)/nameFilter.txt) 
mafft_plus:	$(patsubst %, $(Dir)/mafft_plus/%.fa,	  $(basename $(NAMELIST.1)))

###################################################

#Mkdir
md:make_dir.sh
	bash $<


#Extract all gapless alignment
extra_nog:extract_no_gap.sh 
	bash $< mapped_cds_mafft

#Create "perfect" simulated alignment
$(Dir)/mafft_ali/%.fa:simulate2.R $(Dir)/mafft_nog/%.fa $(Dir)/mafft_ali/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^)

#Remove the gaps 
rm_gap:rem_gaps_ref.sh $(Dir)/mafft_ali 
	bash $< $(word 2,$^) $(Dir)/mafft_rm

#Translate into amino acid and fix the early stop codon issue. 
$(Dir)/aa_seq/%.fa:../dna_to_aa.2.R $(Dir)/mafft_rm/%.fa $(Dir)/aa_seq/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) ${array}

$(Dir)/aa_seq_fix/%.fa:../fix_early_stop.R $(Dir)/aa_seq/%.fa $(Dir)/aa_seq_fix/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) 


#Re-align the sequence via MAFFT
$(Dir)/mafft_reali/%.fa:$(Dir)/aa_seq_fix/%.fa
	$(MAFFT) $< > $@

#Re-map the gaps back to the mafft_rm. 
$(Dir)/mafft_mapped/%.fa: ../mapping_cds2aa.R  $(Dir)/mafft_reali/%.fa $(Dir)/mafft_rm/%.fa $(Dir)/mafft_mapped/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) ${array}

#update the gaps
$(Dir)/updated_gaps/%.fa: ../update_gap.R  $(Dir)/mafft_mapped/%.fa $(Dir)/updated_gaps/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)

#sw-methods on the gaps
$(Dir)/sw_gaps/%.fa: ../pseudo_seq.R  $(Dir)/updated_gaps/%.fa $(Dir)/sw_gaps/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window)

#create a new nameList  
new:../../mafft.nameFilter.sh
	bash $< 

#Add '+' to mafft_In/ and redirect to mafft_plus/
$(Dir)/mafft_plus/%.fa: ../update_gap_plus.R  $(Dir)/sw_gaps/%.fa $(Dir)/mafft_In/%.fa $(Dir)/mafft_plus/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^)

#Generate Displacement table and pie chart of INDELs of different lengths.
mafft_dis:../pos_align_matrix.R $(Dir)/mafft_plus $(Dir)/sw_gaps 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(Dir)/Res/dis.mafft.txt $(Dir)/Res/dis.mafft.pdf $(window) 





clean:
	@rm nogaps.csv
	@rm $(Dir)/mafft_nog/*
	@rm $(Dir)/mafft_ali/*
	@rm $(Dir)/mafft_rm/*
	@rm $(Dir)/mafft_reali/*
	@rm $(Dir)/mafft_mapped/*
	@rm $(Dir)/updated_gaps/*
	@rm $(Dir)/sw_gaps/*
	@rm $(Dir)/mafft_In/*
	@rm $(Dir)/mafft_plus/*

