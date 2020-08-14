default: all
all: 	 Fig
Fig:     mafft_pdf prank_pdf

.PHONY: default all Fig \
	    mafft_pdf prank_pdf

RSCRIPT = Rscript --vanilla 
script  = ../../Script/distance_plot.R

fig.1 = Figure/dis.mafft.pdf
fig.2 = Figure/dis.prank.pdf

#QC --- measure the genetic distance
mafft_pdf: $(script) Alignments/mafft_out mapped_cds_mafft $(INPUT)
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(fig.1) 

prank_pdf: $(script) Alignments/prank_out mapped_cds_prank $(INPUT)
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(fig.2) 


