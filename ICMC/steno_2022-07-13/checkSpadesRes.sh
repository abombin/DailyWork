for i in $(cat newdir.txt); do
if [ -f ./gtdbtk/input/"$i".fasta ]; then
    echo "$i exists" >> ./gtdbtk/completeSpades.txt
    echo "$i exists"
else 
    echo "$i does failed." >> ./gtdbtk/failedSpades.txt
    echo "$i does failed."
fi
done
