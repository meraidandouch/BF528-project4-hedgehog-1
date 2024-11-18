for acc in SRR3879604 SRR3879605 SRR3879606; do
	fasterq-dump $acc & 
done 
wait 
