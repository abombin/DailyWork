# updates from Hflu version
# added numeric check for trimgalore results
# added verbal check for spades results
# use gtdbtk-2.1.0 reduced nessesary instance from m5.16xlarge to m6a.4xlarge (16CPUs 64G mem)
# specified to use 8 cpus for spades

sh ./programs/prepData.sh

cd process_par; ls -d */ | parallel -j 23 'cd {} && sh ../../programs/trimgalore.sh'; cd ../

echo "Finished Trimgalore"

# check number of created files output
ls process_par/*/*_val_1.fq.gz | wc -l; ls process_par/*/*_val_2.fq.gz | wc -l

sleep 10s


cd process_par; ls -d */ | parallel -j 4 'cd {} && sh ../../programs/parSpades.sh'; cd ../

echo "Finished Spades"
sleep 10s

sh ./programs/checkSpadesRes.sh

bash -i ./programs/gtdbtk.sh

echo "Finished GTDB-Tk"

Rscript --vanilla ./programs/findMismatchGtdk.R

bash -i ./programs/forBactopiaGtdbtk.sh

bash -i ./programs/bactopia.sh

echo "Finished Bactopia Run"
sleep 10s

bash -i ./programs/pangenome.sh

sh ./programs/pack.sh

sh ./programs/transferAws.sh









