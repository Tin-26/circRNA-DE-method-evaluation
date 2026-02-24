import pandas as pd
import sys
import glob
from intervaltree import Interval, IntervalTree
from tqdm import tqdm
import os
import re

circRNA_folder = sys.argv[1]

ce2_result_files = glob.glob(os.path.join(circRNA_folder, '*', '*_CE2_quant.txt'))
ce2_fsj_files = glob.glob(os.path.join(circRNA_folder, '*', '*', 'csj_count.txt'))

CE2_columns = ['chrom','start','end','name','score','strand','thickStart','thickEnd',
               'itemRgb','exonCount','exonSizes','exonOffsets','readNumber','circType',
               'geneName','isoformName','index','flankIntron','FPBcirc','FPBlinear','CIRCscore']

csj_columns = ['Chr','Start','End','Strand','Metadata','ReadCount','ReadIDs']

# --- Load BSJ/circRNA table ---
bsj_dfs = []
for file in ce2_result_files:
    df = pd.read_csv(file, sep='\t', header=None, names=CE2_columns, usecols=['chrom','start','end','readNumber','geneName'])
    df['circ_id'] = df['chrom'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    sample_name = os.path.basename(file).split('_')[0]
    df = df[['circ_id','geneName','readNumber']].rename(columns={'readNumber': sample_name})
    bsj_dfs.append(df)

bsj_df = bsj_dfs[0][['circ_id','geneName']].copy()
for df in bsj_dfs:
    bsj_df = pd.merge(bsj_df, df, on=['circ_id','geneName'], how='outer')
bsj_df = bsj_df.fillna(0)
bsj_df.iloc[:,2:] = bsj_df.iloc[:,2:].astype(int)

# --- Parse circRNA coordinates ---
def parse_coords(df, id_col, chrom_col='Chr', start_col='Start', end_col='End'):
    chroms, starts, ends = [], [], []
    for x in df[id_col]:
        chrom, coords = x.split(':')
        start, end = coords.split('-')
        chroms.append(chrom)
        starts.append(int(start))
        ends.append(int(end))
    df[chrom_col] = chroms
    df[start_col] = starts
    df[end_col] = ends
    return df

bsj_df = parse_coords(bsj_df.copy(), 'circ_id')

# --- Process FSJ files and create circRNA x sample matrix ---
fsj_matrix = pd.DataFrame({'circ_id': bsj_df['circ_id'], 'GeneName': bsj_df['geneName']})

for file in tqdm(ce2_fsj_files, desc="Processing FSJ files"):
    sample_name = os.path.basename(os.path.dirname(file))

    df = pd.read_csv(file, sep='\t', header=None, names=csj_columns)
    df['ReadCount'] = pd.to_numeric(df['ReadCount'], errors='coerce').fillna(0).astype(int)
    df['GeneName'] = df['Metadata'].astype(str).str.split(':').str[0]
    df['ID'] = df['Chr'].astype(str)+':'+df['Start'].astype(str)+'-'+df['End'].astype(str)
    df = df[['ID','GeneName','Chr','Start','End','ReadCount']]
    
    # Build interval trees per chromosome
    trees = {}
    for chrom in df['Chr'].unique():
        trees[chrom] = IntervalTree(Interval(row.Start, row.End, row.ReadCount) for row in df[df['Chr']==chrom].itertuples(index=False))
    
    counts = []
    for row in bsj_df.itertuples(index=False):
        chrom = row.Chr
        if chrom not in trees:
            counts.append(0)
            continue
        start_hits = trees[chrom].overlap(row.Start, row.Start+1)
        end_hits = trees[chrom].overlap(row.End, row.End+1)
        total_hits = list(start_hits) + list(end_hits)
        counts.append(sum([h.data for h in total_hits]) if total_hits else 0)
    
    fsj_matrix[sample_name] = counts

# --- Final FSJ matrix (circRNAs x samples) ---
fsj_matrix.fillna(0, inplace=True)
fsj_matrix.iloc[:, 2:] = fsj_matrix.iloc[:, 2:].astype(int)

# --- Save BSJ table ---
bsj_outfile = os.path.join(circRNA_folder, "CE2-BSJ_Matrix.txt")
bsj_df.drop(['Chr', 'Start', 'End'], axis=1, inplace=True)
bsj_df.to_csv(bsj_outfile, index=False, sep='\t')
print(f"BSJ table saved to {bsj_outfile}")

# --- Save FSJ matrix ---
fsj_outfile = os.path.join(circRNA_folder, "CE2-FSJ_Matrix.tx")
fsj_matrix.to_csv(fsj_outfile, index=False, sep='\t')
print(f"FSJ matrix saved to {fsj_outfile}")
