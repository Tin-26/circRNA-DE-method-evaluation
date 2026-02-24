import pandas as pd
import sys
import glob
import os
from tqdm import tqdm

circRNA_folder = sys.argv[1]

ciri_files = glob.glob(os.path.join(circRNA_folder, '*.ciri'))

ciri_candidates = {}
for file in tqdm(ciri_files, total=len(ciri_files), desc="Collecting CIRI3 candidates"):
    file_name = os.path.basename(file).split('_')[1]
    sample_name = file_name.split('.')[0]
    
    df = pd.read_csv(file, sep='\t')
    
    sample_dict = {}
    for _, row in df.iterrows():
        candidate_id = row['circRNA_ID']
        j_reads = row['#junction_reads']
        non_j_reads = row['#non_junction_reads']
        gene_id = row['gene_id']
        alignment_scores = row['SM_MS_SMS']
        ciri_scores = row['Score']
        ratio = row['junction_reads_ratio']
        sample_dict[candidate_id] = [j_reads, non_j_reads, gene_id, alignment_scores, ciri_scores, ratio]
        
    ciri_candidates[sample_name] = sample_dict


unique_ciri_candidates = set()

for sample_dict in ciri_candidates.values():
    unique_ciri_candidates.update(sample_dict.keys())

ciri_gene_ids = {}
for sample_dict in ciri_candidates.values():
    for cid, values in sample_dict.items():
        if cid not in ciri_gene_ids:
            ciri_gene_ids[cid] = values[2]

##CIRI##
bsj_df = pd.DataFrame({'circRNA_ID': sorted(unique_ciri_candidates)})
non_bsj_df = bsj_df.copy()

circ_ids = sorted(unique_ciri_candidates)
# BSJ
bsj_cols = {sample_name: [sample_dict.get(cid, [0,0])[0] for cid in circ_ids] 
            for sample_name, sample_dict in tqdm(ciri_candidates.items(), total=len(ciri_candidates), desc="Creating BSJ matrix")}
bsj_df = pd.DataFrame(bsj_cols, index=circ_ids)
bsj_df.index.name = 'circRNA_ID'
bsj_df.reset_index(inplace=True)
bsj_df['gene_id'] = bsj_df['circRNA_ID'].map(ciri_gene_ids)
# Non-BSJ
non_bsj_cols = {sample_name: [sample_dict.get(cid, [0,0])[1] for cid in circ_ids] 
                for sample_name, sample_dict in tqdm(ciri_candidates.items(), total=len(ciri_candidates), desc="Creating FSJ matrix")}
non_bsj_df = pd.DataFrame(non_bsj_cols, index=circ_ids)
non_bsj_df.index.name = 'circRNA_ID'
non_bsj_df.reset_index(inplace=True)

bsj_df.to_csv(os.path.join(circRNA_folder, 'CIRI-Candidate.BSJ_Matrix'), sep='\t', index=False)
non_bsj_df.to_csv(os.path.join(circRNA_folder, 'CIRI-Candidate.FSJ_Matrix'), sep='\t', index=False)



