#Create a file with the name of the pairwise aligned sequences that DON'T contain gaps

dir=../../Raw_data.2.outgroups


if [ $# -lt 1 ]
then
	echo "At least one argument is required"
fi

for file in ${dir}/$1/*
do
	c=$(grep -c '-' ${file})
	if [ ${c} -eq 0 ]
	then
		echo $(basename ${file}) >> nogaps.csv
	fi
done

########extract the first 2000 files


head -2000 nogaps.csv | while read F
do
	cp $dir/$1/$F $dir/mafft_sim/mafft_nog/
done 
