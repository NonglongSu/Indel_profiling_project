#!/bin/bash

#Create a simulated dna sequence for alignment

default:all

all: md extra_nog sim_align up_gap sw_map new up_plus dis

.PHONY: md extra_nog sim_align up_gap sw_map new up_plus dis 
.PHONY: clean

################################################

RSCRIPT=Rscript --vanilla 
array="c('mouse','rat')"

Dir=../../Raw_data.2.outgroups/mafft_sim
NAMELIST:=$(shell head -2000 nogaps.csv)
window = 6

sim_align:   $(patsubst %, $(Dir)/mafft_ali/%.fa,      $(basename $(NAMELIST)))
up_gap:      $(patsubst %, $(Dir)/updated_gaps/%.fa,   $(basename $(NAMELIST)))
sw_map:      $(patsubst %, $(Dir)/sw_gaps/%.fa,        $(basename $(NAMELIST)))

NAMELIST.1:=$(shell cat $(Dir)/nameFilter.txt) 
up_plus:     $(patsubst %, $(Dir)/mafft_plus/%.fa,     $(basename $(NAMELIST.1)))

###################################################


##Extract all gapless alignment
#extra_nog:extract_no_gap.sh 
#	bash $< mapped_cds_mafft
#
##Create "perfect" simulated alignment
#$(Dir)/mafft_ali/%.fa:simulate2.R $(Dir)/mafft_nog/%.fa $(Dir)/mafft_ali/
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^)

#update the gaps
$(Dir)/updated_gaps/%.fa: ../update_gap.R  $(Dir)/mafft_ali/%.fa $(Dir)/updated_gaps/ 
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
dis:../pos_align_matrix.R $(Dir)/mafft_plus $(Dir)/sw_gaps 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(Dir)/Res/dis.sim.txt $(Dir)/Res/dis.sim.pdf $(window) 





clean:
	@rm nogaps.csv
	@rm $(Dir)/mafft_nog/*
