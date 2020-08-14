default: all
all: 	 Fig 
Fig: 	 mafft_phase  mafft_dis mafft_edge  mafft_edge_report mafft_phase_sep mafft_phase_eff 

.PHONY: default all Fig 
.PHONY: mafft_phase  mafft_dis mafft_edge  mafft_edge_report  mafft_phase_sep mafft_phase_eff

################################################PART I
RSCRIPT = Rscript --vanilla 
sub1    = Figure/mafft_edge

scr1    = ../../Script/indel_phasing.R
scr2    = ../../Script/pos_align_matrix.R
scr3    = ../../Script/edge_bias.R
scr4    = ../../Script/edge_report.R
scr5    = ../../Script/in_del_plot.R 
scr6    = ../../Script/indel_phasing_eff.R

fm1  = Results/phase.mafft.txt
fm2  = Results/dis.mafft.txt
fm3  = Results/phase.mafft.eff.txt

Fm1  = Figure/phase.mafft.pdf
Fm2  = Figure/dis.mafft.pdf
Fm3  = Figure/mafft.edge.report.html
Fm4  = Figure/phase.mafft.sef.pdf
Fm5  = Figure/phase.mafft.eff.pdf


#sub2 = Figure/prank_edge
#fp1  = Results/phase.prank.txt
#fp2  = Results/dis.prank.txt
#fp3  = Results/phae.mafft.eff.txt

#Fp1 = Figure/phase.prank.pdf
#Fp2 = Figure/dis.prank.pdf
#Fp3 = Figure/prank.edge.report.html
#Fp4 = Figure/phase.prank.sep.pdf
#Fp5 = Figure/phase.prank.eff.pdf



#####################################################MAFFT

#Count the proportion of phase 0,1,2 indels
mafft_phase: $(scr1) Mafft/mapped_cds 
	$(RSCRIPT) $< $(word 2,$^) $(fm1) $(Fm1) $(window)

#Generate Displacement table and pie chart of INDELs of different lengths.
mafft_dis:$(scr2) Mafft/updated_cds Mafft/mapped_cds
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(fm2) $(Fm2) $(window) 


#Generate edge_bias table
mafft_edge:$(scr3) Mafft/updated_cds Mafft/mapped_cds Figure/mafft_edge/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(window)

#Generate the edge bias HTML report
mafft_edge_report:$(scr4) $(sub1)/align.ori.l.txt $(sub1)/align.better.l.txt $(sub1)/align.ori.r.txt $(sub1)/align.better.r.txt 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) $(Fm3)


#Count the indel phases of each focal-species
mafft_phase_sep:$(scr5) Mafft/updated_cds Mafft/mapped_cds  
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(Fm4) $(NAME) $(window)

# Count the effective phase distribution of phase 1,2 Indels
mafft_phase_eff:$(scr6) Mafft/mapped_cds
	$(RSCRIPT) $< $(word 2,$^) $(Fm5)






#######################################################PRANK
#prank_phase:$(scr1) Prank/mapped_cds 
#	$(RSCRIPT) $< $(word 2,$^) $(fp1) $(Fp1)

#prank_dis:$(scr2) Prank/mapped_cds_plus Prank/mapped_cds 
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(fp2) $(Fp2) $(window)

#prank_edge:$(scr3) Prank/mapped_cds_plus Prank/mapped_cds Figure/prank_edge/ 
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(window)

#prank_edge_report:$(scr4) $(sub2)/align.ori.l.txt $(sub2)/align.better.l.txt $(sub2)/align.ori.r.txt $(sub2)/align.better.r.txt
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) $(Fp3)

#prank_phase_sep:$(scr5) Prank/updated_cds Prank/mapped_cds
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(Fp4) $(NAME) $(window)

#prank_phase_eff:$(scr6) Prank/mapped_cds
#	$(RSCRIPT) $< $(word 2,$^) $(Fp5)


###################################

#mafft_edge_report:$(scr4)
#        $(RSCRIPT) "rmarkdown::render('$<',output_file='Mafft_edge/Edge.bias.report.html')"
#
#prank_edge_report:$(scr4)
#        $(RSCRIPT) "rmarkdown::render('$<',output_file='Prank_edge/Edge.bias.report.html')"

