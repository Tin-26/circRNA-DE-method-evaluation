import os
import sys
import subprocess
from datetime import datetime
from tqdm import tqdm

project_dir = sys.argv[1]
scripts_dir = sys.argv[2]
ref_fa = sys.argv[3]
ref_gtf = sys.argv[4]
flat_ref = sys.argv[5]
hisat_idx = sys.argv[6]
threads = sys.argv[7]

## Utilities
def run_cmd(cmd):
	try:
		process = subprocess.Popen(
					cmd,
					stdout=subprocess.PIPE,
					stderr=subprocess.STDOUT,
					text=True,
					bufsize=1
			)
		for line in iter(process.stdout.readline, ''):
			sys.stdout.write(line)
			sys.stdout.flush()
		process.stdout.close()
		return_code = process.wait()
		if return_code != 0:
			print(f"Command failed with return code: {return_code}")
	except Exception as e:
		print(f"Error: {e}")

def log(msg):
		print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}")

### Scripts
fastp_script = os.path.join(scripts_dir, "fastp.py")
fastqc_script = os.path.join(scripts_dir, "fastqc.py")
bwa_map_script = os.path.join(scripts_dir, "bwa_map.py") # Mapping for CIRI3 & CE2
ciri_script = os.path.join(scripts_dir, "CIRI_ID.py") # ID with CIRI
ce2_script = os.path.join(scripts_dir, "CE2_ID.py") # ID with CE2
matrix_script = os.path.join(scripts_dir, "CIRI_matrix.py") # Creates BSJ,FSJ matrices from CIRI results
ce2_matrix_script = os.path.join(scripts_dir, "CE2_matrix.py") # BSJ, FSJ matrix from CE2

### Folders
raw_fastq_dir = os.path.join(project_dir, "raw")
cleaned_fastq_dir = os.path.join(project_dir, "cleaned")
fastp_reports_dir = os.path.join(project_dir, "FASTP_reports")
bwa_output_folder = os.path.join(project_dir, "aligned")
ciri_output_folder = os.path.join(project_dir, "circRNA_outs", "CIRI3")
ce2_base_folder = os.path.join(project_dir, "circRNA_outs", "CE2")
os.makedirs(ciri_output_folder, exist_ok=True)
os.makedirs(bwa_output_folder, exist_ok=True)

## Step 1. ##
fastp_cmd = ['python3', fastp_script, raw_fastq_dir, cleaned_fastq_dir, fastp_reports_dir,'no']
log("## Cleaning fastq files ##")
run_cmd(fastp_cmd)

## Step 2. ##
fastqc_cmd = ['python3', fastqc_script, project_dir, threads]
log("## Running QC ##")
run_cmd(fastqc_cmd)

file_list = sorted(set(os.listdir(cleaned_fastq_dir)))

with tqdm(total = len(file_list)/2, desc = "Processing samples", dynamic_ncols=True, leave=False) as pbar: 
	for sample in file_list:
		name = sample.split('_')[0]
		forward_reads = os.path.join(cleaned_fastq_dir, sample)
		reverse_reads = os.path.join(cleaned_fastq_dir, sample.replace("_R1_cleaned", "_R2_cleaned"))
	
		## Step 3. ##
	
		sam_file = os.path.join(bwa_output_folder, f"{name}.sam")
		viewed_bam = os.path.join(bwa_output_folder, f"{name}_viewed.bam")
		sorted_bam = os.path.join(bwa_output_folder, f"{name}.bam")
	
		hisat_bam = os.path.join(bwa_output_folder, f"HT2-{name}.bam")
		hisat_cmd = (
				f"hisat2 -x {hisat_idx} -1 {forward_reads} -2 {reverse_reads} -p {threads} "
	    		f"| samtools view -@ {threads} -b -F 4 "
	    		f"| samtools sort -@ {threads} -o {hisat_bam}"
			)
	
		align_cmd = [
			"bwa", "mem",
			"-t", threads,
			"-T", "19",
			ref_fa,
			forward_reads,
			reverse_reads
		]
	
		cmd_view = ["samtools", "view", "-@", threads, "-S", "-b", sam_file, "-o", viewed_bam]
		cmd_sort = ["samtools", "sort", "-@", threads, viewed_bam, "-o", sorted_bam]
		cmd_index = ["samtools", "index", "-@", threads, sorted_bam]
		hisat_index_cmd = ["samtools", "index", "-@", threads, hisat_bam]
		if not os.path.exists(sam_file):
			print(f"[{name}] Aligning...")
			with open(sam_file, "w") as f:
				process = subprocess.Popen(align_cmd, stdout=f)
				process.communicate()
	
			print(f"[{name}] Converting SAM to BAM...")
		if not os.path.exists(sorted_bam):
			run_cmd(cmd_view)
		
			#print(f"[{name}] Sorting BAM...")
			run_cmd(cmd_sort)
			#print(f"[{name}] Indexing BAM...")
			run_cmd(cmd_index)

		# Clean up intermediate files
		if os.path.exists(viewed_bam):
			os.remove(viewed_bam)
	#
		if not os.path.exists(hisat_bam):
			try:
				log(f"[{name}] Hisat2 aligning...")
				subprocess.run(hisat_cmd, shell=True, check=True)
			except subprocess.CalledProcessError as e:
				log(f"[ERROR] Something happened during running pipeline: {e}")
		if not os.path.exists(os.path.join(bwa_output_folder, f"{hisat_bam}.bai")):
			run_cmd(hisat_index_cmd)
	
		## Step 4. ##	
		log(f"## Identifying circRNAs with CIRI3 for {name} ##")
		ciri_out = os.path.join(ciri_output_folder, f'CIRI3_{name}.ciri')
		ciri_id = ['java', '-jar', '/media/meteor/FatDawg/Benchmark_Paper/Scripts/CIRI3_Java_18.0.1.jar',
					'-I', sam_file,
					'-O', ciri_out,
					'-F', ref_fa,
					'-A', ref_gtf,
					'-T', threads 
					]
		if not os.path.exists(ciri_out):
			run_cmd(ciri_id)
	
		log(f"## Identifying circRNAs with CE2 for {name} ##")
		ce2_output_folder = os.path.join(project_dir, "circRNA_outs", "CE2", name)
		os.makedirs(ce2_output_folder, exist_ok=True)
		parser_out = f"back_spliced_junction.bed"
		annotate_out = f"{name}_circRNAs.txt"
		container = "/media/meteor/FatDawg/Benchmark_Paper/Scripts/CIRCexplorer2.sif"
		ce2_parse_cmd = ["singularity", "exec",
					"--bind", f"{bwa_output_folder}:/input",
					"--bind", f"{ce2_output_folder}:/out",
					container,
					"bash", "-c", f"cd /out && CIRCexplorer2 parse -t BWA /input/{name}.bam" 
				]
		ce2_annotate_cmd = ["singularity", "exec",
					"--bind", f"{ce2_output_folder}:/app",
					"--bind", f"{os.path.dirname(os.path.abspath(ref_fa))}:/references",
					#"--bind", f"{annotate_out_folder}:/circexplorer2_outs",
					container,
					"CIRCexplorer2", "annotate",
					"-r", "/references/Hsapiens_GRCh38_112_GenePred_reference_flat.txt",
					"-g", "/references/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
					"-b", f"/app/{parser_out}",
					"-o", f"/app/{annotate_out}"
				]
	
		if not os.path.exists(os.path.join(ce2_output_folder, parser_out)):
			log("## [CE2] Parsing results ##")
			run_cmd(ce2_parse_cmd)
		if not os.path.exists(os.path.join(ce2_output_folder, annotate_out)):
			log("## [CE2] Annotating results ##")
			run_cmd(ce2_annotate_cmd)
		
		# CLEAR/CircExplorer3 pipeline
		ce2_quant_out = os.path.join(ce2_output_folder, f"{name}_CE2_quant.txt")
		ce2_quant = ['circ_quant', '-c', os.path.join(ce2_output_folder, annotate_out), '-b', hisat_bam, '-t', '-r', flat_ref, '-o', ce2_quant_out]
		if not os.path.exists(ce2_quant_out):
			log("## [CE2] Quantifying ##")
			run_cmd(ce2_quant)
		#if os.path.exists(sam_file):
		#	os.remove(sam_file)
		pbar.update(1)

## Step 6. ## Create circRNA matrices
create_ciri_matrices = ['python3', matrix_script, ciri_output_folder]
create_ce2_matrices = ['python3', ce2_matrix_script, ce2_base_folder]
ciri_bsj_matrix = os.path.join(ciri_output_folder, 'CIRI-Candidate.BSJ_Matrix')
ciri_fsj_matrix = os.path.join(ciri_output_folder, 'CIRI-Candidate.FSJ_Matrix')
ce2_bsj_matrix = os.path.join(ciri_output_folder, 'CE2-BSJ_Matrix.txt')
ce2_fsj_matrix = os.path.join(ciri_output_folder, 'CE2-FSJ_Matrix.txt')

if not os.path.exists(ciri_bsj_matrix):
	log(f"## Creating matrices from CIRI3 results ##")
	run_cmd(create_ciri_matrices)
if not os.path.exists(ce2_bsj_matrix):
	log(f"## Creating matrices from CE2 results ##")
	run_cmd(create_ce2_matrices)

log(f"## All samples are done!")
