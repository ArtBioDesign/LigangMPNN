#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName     :run.py
# @Time         :2024/04/18 19:47:21
# @Author       :YangChunhe
# @Email        :2393492851@qq.com
# @Description  :file content
import os
import sys
import pandas as pd
import time
import json
import argparse
  
# warnings.filterwarnings('ignore')

def pdb_to_csv(site_pdb_path,outpudir):
    df = get_rn_from_pdb(site_pdb_path, df_name="ATOM", columns_name=["chain_id", "residue_number"])
    df = df.rename(columns={"chain_id": "chain", "residue_number": "pos"})
    name = os.path.splitext(os.path.basename(site_pdb_path))[0]
    site_file = os.path.join(outpudir, f"{name}.csv")
    df.to_csv(site_file, index=False)
    return site_file


def get_rn_from_pdb(pdb_path, df_name="ATOM", columns_name=["chain_id", "residue_number"]):
    from biopandas.pdb import PandasPdb
    ppdb = PandasPdb().read_pdb(pdb_path)
    df = ppdb.df["ATOM"][columns_name].drop_duplicates()
    return df
    

def write_num(out_dir,time,num):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    filename = os.path.join(out_dir ,"time_num.json")
    data={
        "time":time,
        "num":num
    }
    try:
        with open(filename, 'w') as f:
            json.dump(data, f, indent=4)
        print(f"Data has been written to {filename}")
    except Exception as e:
        print(f"An error occurred while writing to {filename}: {e}")
    
def run(site_pdb_path, out_dir):
    max_workers = 4
    s0 = time.time()
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    site_csv_file = pdb_to_csv(site_pdb_path, out_dir)
    task1_df = pd.read_csv(site_csv_file)
    task1_df['pos'] = task1_df['pos'].astype(str)
    task1_df["chain_pos"] = task1_df['chain'] + task1_df['pos']
    residues = list(set(task1_df["chain_pos"]))
    s1 = time.time()
    print("compute time cost: %f s"%(s1-s0))
    write_num(out_dir,s1-s0,len(residues))
        
def main():
    parser = argparse.ArgumentParser(description='check site counts')
    parser.add_argument('--site_pdb', default="./inputs/8A-mix43.pdb", help='site pdb file')
    parser.add_argument('--out_dir', default='./outputs/BsYjiC/', help='output dir')
    args = parser.parse_args()
    run( args.site_pdb, args.out_dir )

   

if __name__ == "__main__":
    main()
