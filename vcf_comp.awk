function rm_genotype(s) 
{
	return gensub("=.*$","",s,s)
}
function rm_samples(s)
{
	split(s, r, "=")
	if (r[2] != "./.") {
	return r[2]
	} else {
		return "NA"
	}
}

NR==1{
	min = split($0,samples,"\t");

	printf "CHR" "\t" "POS" "\t" "TYPE" "\t" "REF" "\t" "ALT" sample_list; \
	for ( i=6; i<=NF; i++ ) {
		printf("\t%s",rm_genotype(samples[i]));
	}
	printf "\n"
}

{
	printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 
	for ( i=6; i<=NF; i++ ) {
	printf("\t%s",rm_samples($i));
	}
	printf "\n"
}