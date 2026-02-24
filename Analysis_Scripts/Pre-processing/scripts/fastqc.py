import os
import subprocess
import glob
import sys

project_folder = sys.argv[1]
threads = sys.argv[2]

raw_folder = os.path.join(project_folder, 'raw')
cleaned_folder = os.path.join(project_folder, 'cleaned')
output_folder = os.path.join(project_folder, 'FASTQC_reports')

pre_trim_folder = os.path.join(output_folder, 'pre_trim')
post_trim_folder = os.path.join(output_folder, 'post_trim')

os.makedirs(pre_trim_folder, exist_ok=True)
os.makedirs(post_trim_folder, exist_ok=True)

def has_fastqc_output(fastq_file, output_folder):
    base_name = os.path.basename(fastq_file)
    if base_name.endswith('.gz'):
        base_name = base_name[:-3]
    if base_name.endswith('.fastq') or base_name.endswith('.fq'):
        base_name = base_name.rsplit('.', 1)[0]

    expected_html = os.path.join(output_folder, base_name + '_fastqc.html')
    expected_zip = os.path.join(output_folder, base_name + '_fastqc.zip')
    return os.path.exists(expected_html) and os.path.exists(expected_zip)

def run_fastqc(input_folder, output_folder, threads=28):
    fastq_files = glob.glob(os.path.join(input_folder, '*.fastq')) + \
                  glob.glob(os.path.join(input_folder, '*.fq')) + \
                  glob.glob(os.path.join(input_folder, '*.fastq.gz')) + \
                  glob.glob(os.path.join(input_folder, '*.fq.gz'))

    if not fastq_files:
        print(f"No FASTQ files found in {input_folder}")
        return

    to_process = [f for f in fastq_files if not has_fastqc_output(f, output_folder)]

    if not to_process:
        print(f"All FASTQ files in {input_folder} already have FastQC reports. Skipping.")
        return

    print(f"Running FastQC on {len(to_process)} files in {input_folder}...")

    cmd = [
        "fastqc",
        *to_process,
        "--outdir", output_folder,
        "--threads", str(threads)
    ]

    subprocess.run(cmd, check=True)

print("Running FastQC on raw files...")
run_fastqc(raw_folder, pre_trim_folder, threads)
multiqc_cmd_pre = ['multiqc', pre_trim_folder, '--outdir', output_folder, '--filename', 'Pre_trim_QC']
if not os.path.exists(os.path.join(output_folder, "Pre_trim_QC.html")):
    subprocess.run(multiqc_cmd_pre, check=True)
print("Running FastQC on cleaned files...")
run_fastqc(cleaned_folder, post_trim_folder, threads)
multiqc_cmd_post = ['multiqc', post_trim_folder, '--outdir', output_folder, '--filename', 'Post_trim_QC']
if not os.path.exists(os.path.join(output_folder, "Post_trim_QC.html")):
    subprocess.run(multiqc_cmd_post, check=True)
print("All FastQC analyses completed.")
