#remove all gaps from previous "perfect" alignment

FROM=$1
TO=$2


mkdir -p ${TO}

for file in ${FROM}/*
do
	cat ${file} | tr -d '-' >  ${TO}/$(basename ${file})
done
