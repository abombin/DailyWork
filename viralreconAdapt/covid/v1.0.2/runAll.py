import subprocess

# download files highly situational so commented out
#subprocess.run("python programs/download.py", shell=True)

# make input .csv file
subprocess.run("python programs/inputFile.py", shell=True)

# run Viralrecon
subprocess.run("python programs/runViralrecon.py", shell=True)

# Store results in S3 bucket
subprocess.run("python programs/results2S3.py", shell=True)