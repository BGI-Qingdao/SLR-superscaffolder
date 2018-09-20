
mkdir output
mkdir logs
mkdir tmps

mv xxx.scaff_seqs output/ || echo xxx.scaff_seqs is not exist!!!

mv log_* logs/
mv xxx.* tmps/

echo "info : output scaff is in output/scaff_seqs "
echo "info : logs are all in logs/*"
echo "info : you can delete folder tmps by : rm -rf tmps "
