#!/usr/bin/env python3

import os
import subprocess
import sys

def run_fastp(input1, input2, output1, output2, report_html, report_json):
    command = [
        "fastp",
        "-i", input1,
        "-I", input2,
        "-o", output1,
        "-O", output2,
        "-w", "16", # number of threads used
        "--qualified_quality_phred", "28", # Sets the Phred score threshold for a base to be considered “qualified”! Phred 28 ≈ 0.16% error rate -> High quality reads
        "--cut_front", "3", # Trim 3 bases from the 5' end of each read. ->
        "--cut_tail", "3", # Trim 3 bases from the 3' end of each read.
        "--cut_right", "4", # Enables sliding window cut, trim when window fails quality criteria. (similar to trimmomatic) 
        "--cut_window_size", "4", # window size is 4 bases
        "--cut_mean_quality", "15", # mean quality of the 4 bases must be <= 15
        "--length_required", "36", # Minimum length of the read
        "--low_complexity_filter", #Does the filtering of low-complexity reads (e.g. homopolymers, repeats).
        "-x", # poly X tail trimming
        "-G", # poly G tail trimming     with both enabled G trim goes first, then everything else
        "--low_complexity_filter", # complexity filter:  the percentage of base that is different from its next base
        "--complexity_threshold", "30", # 30% complexity is required
        "--detect_adapter_for_pe", # automatic adapter trimming
        "--html", report_html, # report options
        "--json", report_json
    ]
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] fastp failed: {e}")

def run_SE_fastp(input1, output1, report_html, report_json):
    se_command = [
        "fastp",
        "-i", input1,
        "-o", output1,
        "-w", "16",
        "--qualified_quality_phred", "28",
        "--cut_front", "3",
        "--cut_tail", "3",
        "--cut_right", "4",
        "--cut_window_size", "4",
        "--cut_mean_quality", "15",
        "--length_required", "36",
        "--low_complexity_filter",
        "--complexity_threshold", "30",
        "--html", report_html,
        "--json", report_json
    ]
    try:
        subprocess.run(se_command, check=True, close_fds=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] fastp failed: {e}")

def process_fastq_files(input_folder, output_folder, report_folder, is_se):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if not os.path.exists(report_folder):
        os.makedirs(report_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith("_R1.fastq.gz") and os.path.isfile(os.path.join(input_folder, filename)):
            sample_name = filename.replace("_R1.fastq.gz", "")
            input1_path = os.path.join(input_folder, filename)

            if is_se == 'no':
                paired_filename = filename.replace("_R1.fastq.gz", "_R2.fastq.gz")
                input2_path = os.path.join(input_folder, paired_filename)

                if os.path.isfile(input2_path):
                    output1 = os.path.join(output_folder, filename.replace("_R1.fastq.gz", "_R1_cleaned.fastq.gz"))
                    output2 = os.path.join(output_folder, paired_filename.replace("_R2.fastq.gz", "_R2_cleaned.fastq.gz"))
                    report_html = os.path.join(report_folder, filename.replace("_R1.fastq.gz", "_fastp_report.html"))
                    report_json = os.path.join(report_folder, filename.replace("_R1.fastq.gz", "_fastp_report.json"))

                    if not os.path.exists(report_html):
                        print(f"[INFO] Paired-end QC in progress for {sample_name}")
                        run_fastp(input1_path, input2_path, output1, output2, report_html, report_json)
                        print(f"[DONE] QC completed for {sample_name}, report: {report_html}")
                else:
                    print(f"[WARN] Missing pair file for {filename}: expected {paired_filename}")
            else:  # Single-end
                output = os.path.join(output_folder, filename.replace("_R1.fastq.gz", "_R1_cleaned.fastq.gz"))
                report_html = os.path.join(report_folder, filename.replace("_R1.fastq.gz", "_fastp_report.html"))
                report_json = os.path.join(report_folder, filename.replace("_R1.fastq.gz", "_fastp_report.json"))

                if not os.path.exists(report_html):
                    print(f"[INFO] Single-end QC in progress for {sample_name}")
                    run_SE_fastp(input1_path, output, report_html, report_json)
                    print(f"[DONE] QC completed for {sample_name}, report: {report_html}")

if __name__ == "__main__":
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    report_folder = sys.argv[3]
    is_se = sys.argv[4].strip().lower()

    process_fastq_files(input_folder, output_folder, report_folder, is_se)
