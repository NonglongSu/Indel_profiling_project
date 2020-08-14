
cat gencode.v31.annotation.gtf | grep -v '^#' | grep -E '(cds_start_NF|cds_end_NF)' | awk '$3=="CDS"' | cut -f9 \
	| awk '{print $4}' | cut -d '"' -f2 | cut -d. -f1 | sort | uniq -c | awk '{print $2}' > trans_ID_NF.txt
