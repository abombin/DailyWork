bash -i ./programs/renameFiles.sh

bash -i ./programs/bactopia.sh # needs setting of bacteria name

bash -i ./programs/mlst.sh

Rscript --vanilla ./programs/sumMlst.R

sh ./programs/results2S3.sh