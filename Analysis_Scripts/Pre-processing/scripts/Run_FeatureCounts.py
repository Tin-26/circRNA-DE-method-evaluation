import os
import sys
import subprocess
from datetime import datetime
import glob
import pandas as pd

gtf = sys.argv[1]
input_dir = sys.argv[2]
output_dir = sys.argv[3]
threads = sys.argv[4]
today = datetime.today().strftime('%Y-%m-%d')
output_file = os.path.join(output_dir, f'featureCounts_matrix.txt')
bam_files = glob.glob(os.path.join(input_dir, '*/*.bam'))
renamed_bam_files = []

def run_cmd(cmd):
    """
    Run a command as a subprocess and print its output in real-time.

    :param cmd: List of command arguments to execute.
    :return: None
    """
    try:
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True, 
            bufsize=1,
        )

        for line in process.stdout:
            print(line.strip())

        for line in process.stderr:
            print(line.strip())

        process.wait()

        if process.returncode != 0:
            print(f"Command failed with return code: {process.returncode}")

    except Exception as e:
        print(f"An error occurred: {e}")

for bam in bam_files:
    if "Aligned.sortedByCoords.out" in bam:
        new_bam = bam.replace("Aligned.sortedByCoords.out", "")
        os.rename(bam, new_bam)  # Rename the file
        renamed_bam_files.append(new_bam)
    else:
        renamed_bam_files.append(bam)

if not bam_files:
    print("No BAM files found in the specified directory.")
    sys.exit(1)

cmd = ['featureCounts',
       '-p',
       '--countReadPairs',
       '-T', str(threads),
       '-a', gtf,
       '-o', output_file] + renamed_bam_files
       
run_cmd(cmd)

print("Processing output...")
with open(output_file, 'r') as txt_file:
    lines = txt_file.readlines() 

lines = [line for line in lines if not line.startswith("#")]

with open(output_file, 'w') as txt_file:
    txt_file.writelines(lines)

csv = pd.read_csv(output_file, sep='\t', header = [0], index_col=[0]) 

new_columns = [os.path.basename(col).replace("Aligned.sortedByCoord.out.bam", "") for col in csv.columns]

csv.columns = new_columns
csv.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], inplace=True, axis=1)
csv.to_csv(output_file, sep='\t')
print("All done!")
