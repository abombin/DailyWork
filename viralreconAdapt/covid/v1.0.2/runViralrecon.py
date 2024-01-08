import subprocess
import glob

inputList = glob.glob("*.csv")

def runViralrecon(inputList, trim, protocol):
    for inputFile in inputList:
        runtype = inputFile.replace("_input.csv", "")
        if protocol == "amplicon":
            cmd_str = f"/home/ubuntu/nextflow run nf-core/viralrecon -r 2.6.0 \
                --max_cpus 94 \
                --max_memory '370.GB' \
                --input {inputFile} \
                --outdir {runtype}_output_trim_{trim} \
                --platform illumina \
                --protocol amplicon \
                --genome 'MN908947.3' \
                --kraken2_db ../../references/kraken2-human-db \
                --skip_kraken2 \
                --save_reference false \
                --variant_caller 'ivar' \
                --primer_bed ../../references/swift_refv3_primers.bed \
                --primer_left_suffix '_LEFT' \
                --primer_right_suffix '_RIGHT' \
                --ivar_trim_offset {trim} \
                --skip_assembly \
                -c /home/ubuntu/extraVol/Viralrecon/covid/custom.config \
                -profile docker \
                -with-docker nfcore/virarecon"
        subprocess.run(cmd_str, shell=True)

# run function
runViralrecon(inputList, trim=0, protocol="amplicon")

# make a function abvoe modular via splitting command string to substrins and then use " ".join([])
