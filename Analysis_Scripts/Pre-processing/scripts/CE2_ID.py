import os
import subprocess
import sys
import gc

input_dir = sys.argv[1]
reference = sys.argv[2]
gtf = sys.argv[3]
threads = str(sys.argv[4])
container = "/media/meteor/FatDawg/Benchmark_Paper/Scripts/CIRCexplorer2.sif"
def run_cmd(cmd):
	"""
	Run a command as a subprocess and print its output in real-time.
	:param cmd: List of command arguments to execute.
	"""
	try:
		process = subprocess.Popen(
			cmd,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE,
			text=True,
			bufsize=1
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

dir_list = sorted(
    [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]
)


# Process each sample
for sample in dir_list:
	print(f"--- Processing sample: {sample} ---\n")
	name = sample
	input_bam = os.path.join(input_dir, name, f"{name}.bam")
	sam_file = os.path.join(input_dir, name, f"{name}.sam")
	viewed_bam = os.path.join(input_dir, name, f"{name}_viewed.bam")
	output = os.path.join(input_dir, name, "CE2_outs")
	parser_out_folder = os.path.join(output, "Parser_output")
	annotate_out_folder = os.path.join(output, "Annotate_output")
	parser_out = f"back_spliced_junction.bed"
	annotate_out = f"{name}_circRNAs.txt"
	os.makedirs(output, exist_ok = True)

	cmd_view = ["samtools", "view", "-@", threads, "-S", "-b", sam_file, "-o", viewed_bam]
	cmd_sort = ["samtools", "sort", "-@", threads, viewed_bam, "-o", input_bam]
	cmd_index = ["samtools", "index", "-@", threads, input_bam]
	ce2_parse_cmd = ["singularity", "exec",
				"--bind", f"{os.path.join(input_dir, name)}:/input",
				"--bind", f"{output}:/out",
				container,
				"bash", "-c", f"cd /out && CIRCexplorer2 parse -t BWA /input/{name}.bam" 
			]
	ce2_annotate_cmd = ["singularity", "exec",
				"--bind", f"{output}:/app",
				"--bind", f"{os.path.dirname(os.path.abspath(reference))}:/references",
				#"--bind", f"{annotate_out_folder}:/circexplorer2_outs",
				container,
				"CIRCexplorer2", "annotate",
				"-r", "/references/Hsapiens_GRCh38_112_GenePred_reference_flat.txt",
				"-g", "/references/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
				"-b", f"/app/{parser_out}",
				"-o", f"/app/{annotate_out}"
			]

	if not os.path.exists(input_bam):
		run_cmd(cmd_view)
		run_cmd(cmd_sort)
		run_cmd(cmd_index)
	
	if not os.path.exists(os.path.join(output, parser_out)):
		run_cmd(ce2_parse_cmd)
	if not os.path.exists(os.path.join(output, annotate_out)):
		run_cmd(ce2_annotate_cmd)

	# Free up memory
	del ce2_parse_cmd, ce2_annotate_cmd
	if os.path.exists(viewed_bam):
		os.remove(viewed_bam)
	gc.collect()

	print(f"--- Finished sample: {name} ---\n")
