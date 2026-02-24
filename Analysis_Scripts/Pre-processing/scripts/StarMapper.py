#!/bin/python

import os
import sys
import subprocess

genome_path = sys.argv[1]
file_list = os.listdir(sys.argv[2])
file_list.sort()
destination_folder = sys.argv[3]
threads = sys.argv[4]

i = 0
print("============================\nLoading genome...\n============================")
subprocess.run(['STAR', '--genomeLoad', 'LoadAndExit', '--genomeDir', genome_path])
print("============================\nGenome is loaded!\n============================")
for i in  range(0,len(file_list),2):
	fastq_path1 = os.path.join(sys.argv[2], file_list[i])
	fastq_path2 = fastq_path1.replace('_R1_', '_R2_')
	name = file_list[i][:-9]
	name = name.split("_")[0]
	
	output_folder = os.path.join(destination_folder, name)
	if not os.path.exists(output_folder):
		os.makedirs(output_folder, exist_ok=True)
		cmd = ['STAR', '--genomeLoad', 'LoadAndKeep','--runThreadN', str(threads),
		'--genomeDir', genome_path,
	 	'--readFilesIn', fastq_path1, fastq_path2,
	 	'--readFilesCommand', 'zcat',
	 	'--outSAMtype', 'BAM', 'SortedByCoordinate',
	 	'--limitBAMsortRAM', '10000000000',
	 	'--outFileNamePrefix', os.path.join(output_folder, name)
	 	#'--outFilterMatchNmin 18',
	 	
	 	]	
		print('Running:')
		subprocess.run(cmd, check=True)
print("============================\nRemoving loaded genome...\n============================")
subprocess.run(['STAR', '--genomeLoad', 'Remove', '--genomeDir', genome_path])
print("============================\nAll done!\n============================")
