import pandas as pd
import glob 
import os
if __name__=='__main__':

    data_dir = "/Users/amohamed/readdy/Helmut_data/EFs_WT_dark"
    outdir = "/Users/amohamed/readdy/Helmut_data/EFs_WT_dark_proc"
    list_of_files = sorted(glob.glob(f"{data_dir}/*.csv"))

    os.makedirs(outdir, exist_ok=True)
    for csv in list_of_files:
        df = pd.read_csv(csv)
        if not "Label" in df.columns: continue
        filtered_df = df[df['Label'].str.contains('particleC2')]
        filtered_df = filtered_df[['XM', 'YM']]
        filtered_df = filtered_df.rename(columns={'XM': 'X', 'YM': 'Y'})
        filtered_df.to_csv(f"{outdir}/{csv.split('/')[-1].replace('.csv','_PSIIs_.csv')}")


