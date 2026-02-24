import os
import subprocess
import sys
import gc

input_dir = sys.argv[1]
output_dir = sys.argv[2]
reference = sys.argv[3]
gtf = sys.argv[4]
threads = str(sys.argv[5])

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

dir_list = os.listdir(input_dir)


# Process each sample
for sample in dir_list:
	print(f"--- Processing sample: {sample} ---\n")
	name = sample.split(".")[0]
	input_sam = os.path.join(input_dir, f"{name}.sam")
	output = os.path.join(output_dir, name, f"CIRI3_{name}.ciri")
	if not os.path.exists(os.path.join(output_dir, name)):
		os.makedirs(os.path.join(output_dir, name))
	ciri_cmd = ['java', '-jar', '/media/meteor/FatDawg/Benchmark_Paper/Scripts/CIRI3_Java_18.0.1.jar',
				'-I', input_sam,
				'-O', output,
				'-F', reference,
				'-A', gtf,
				'-T', threads 
				]

	if not os.path.exists(output):
		run_cmd(ciri_cmd)

	# Free up memory
	del ciri_cmd
	gc.collect()

	print(f"--- Finished sample: {name} ---\n")
