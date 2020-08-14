#default: all
#all: 	 Fig
Fig: 	mafft_phase prank_phase mafft_dis prank_dis mafft_edge prank_edge mafft_edge_report prank_edge_report

.PHONY: default all Fig
.PHONY: mafft_phase prank_phase mafft_dis prank_dis mafft_edge prank_edge mafft_edge_report prank_edge_report

RSCRIPT=Rscript --vanilla 
script.1=../../../Script/indel_phasing.R
script.2=../../../Script/pos_align_matrix.R
script.3=../../../Script/edge.bias.R
script.4=../../../Script/edge.report.R


window=6

file.1=../Data/Results/phase.mafft.txt
file.2=../Data/Results/phase.prank.txt
file.3=../Data/Results/displacement.mafft.table.txt
file.4=../Data/Results/displacement.prank.table.txt

file.5=mafft_edge/align.ori.l.txt
file.6=mafft_edge/align.better.l.txt
file.7=mafft_edge/align.ori.r.txt
file.8=mafft_edge/align.better.r.txt

file.9 =prank_edge/align.ori.l.txt
file.10=prank_edge/align.better.l.txt
file.11=prank_edge/align.ori.r.txt
file.12=prank_edge/align.better.r.txt

fig.1=phase.mafft.chart.pdf
fig.2=phase.prank.chart.pdf
fig.3=displacement.mafft.pdf
fig.4=displacement.prank.pdf

fig.5=mafft.edge.report.html
fig.6=prank.edge.report.html











#Count the Freq distribution of phase 0,1,2 Indels
mafft_phase:$(script.1) ../Data/Mafft/mapped_cds_Out 
	$(RSCRIPT) $< $(word 2,$^) $(file.1) $(fig.1)

prank_phase:$(script.1) ../Data/Prank/mapped_cds_Out 
	$(RSCRIPT) $< $(word 2,$^) $(file.2) $(fig.2)

#Generate Displacement table and pie chart of INDELs of different lengths.
mafft_dis:$(script.2) ../Data/Mafft/mapped_cds_plus ../Data/Mafft/mapped_cds_Out 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(file.3) $(fig.3) $(window)

prank_dis:$(script.2) ../Data/Prank/mapped_cds_plus ../Data/Prank/mapped_cds_Out 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(file.4) $(fig.4) 


#Generate edge bias table
mafft_edge:$(script.3) ../Data/Mafft/mapped_cds_plus ../Data/Mafft/mapped_cds_Out mafft_edge/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(window)

prank_edge:$(script.3) ../Data/Prank/mapped_cds_plus ../Data/Prank/mapped_cds_Out prank_edge/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(window)


#Generate the edge bias HTML report
mafft_edge_report:$(script.4)
	$(RSCRIPT) $< $(file.5) $(file.6) $(file.7) $(file.8) $(fig.5)

prank_edge_report:$(script.4)
	$(RSCRIPT) $< $(file.9) $(file.10) $(file.11) $(file.12) $(fig.6)






#mafft_edge_report:$(script.4)
#        $(RSCRIPT) "rmarkdown::render('$<',output_file='Mafft_edge/Edge.bias.report.html')"
#
#prank_edge_report:$(script.4)
#        $(RSCRIPT) "rmarkdown::render('$<',output_file='Prank_edge/Edge.bias.report.html')"

