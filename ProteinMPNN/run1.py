#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName     :run.py
# @Time         :2024/04/18 19:47:21
# @Author       :YangChunhe
# @Email        :2393492851@qq.com
# @Description  :file content
import os
import sys
import subprocess
from Bio import SeqIO
import pandas as pd
from pathlib import Path
import multiprocessing
import warnings
import shutil
import gc
import time
import json
import argparse

warnings.filterwarnings('ignore')

def has_files_in_current_directory_and_subdirectories(current_directory):
    directory = Path(current_directory)
    for item in directory.rglob('*'):
        if item.is_file():
            return True
    return False

def get_rn_from_pdb(pdb_path, df_name="ATOM", columns_name=["chain_id", "residue_number"]):
    from biopandas.pdb import PandasPdb
    ppdb = PandasPdb().read_pdb(pdb_path)
    df = ppdb.df["ATOM"][columns_name].drop_duplicates()
    return df

def pdb_to_csv(site_pdb_path,outpudir):
    df = get_rn_from_pdb(site_pdb_path, df_name="ATOM", columns_name=["chain_id", "residue_number"])
    df = df.rename(columns={"chain_id": "chain", "residue_number": "pos"})
    name = os.path.splitext(os.path.basename(site_pdb_path))[0]
    site_file = os.path.join(outpudir, f"{name}.csv")
    df.to_csv(site_file, index=False)
    return site_file

def extract_mutation_info(wildtype, design):
    mutations = [f"{wildtype[i]}{i+1}{design[i]}" for i in range(len(wildtype)) if wildtype[i] != design[i]]
    return mutations[0] if mutations else ""

def process_file(input_file, name):
    data = []
    wildtype = None
    overall_confidence = None
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith(f">{name}, T=0.1"):
                wildtype = next(infile).strip()
            elif line.startswith(f">{name}, id="):
                parts = line.split(", ")
                seq_id = parts[1]
                overall_confidence = float(parts[4].split("=")[1])
                design = next(infile).strip()
                mutation_line = extract_mutation_info(wildtype, design)
                if mutation_line.strip():
                    data.append([seq_id, overall_confidence, mutation_line.strip()])
    return pd.DataFrame(data, columns=['ID', 'Overall Confidence', 'Mutation Info'])

def process_folder(folder_path, output_file, name):
    all_data = []
    for subdir, dirs, files in os.walk(folder_path):
        for dir_name in dirs:
            if dir_name.startswith('redesign') and has_files_in_current_directory_and_subdirectories(os.path.join(subdir, dir_name)):
                seqs_folder = os.path.join(subdir, dir_name, 'seqs')
                input_file_path = os.path.join(seqs_folder, f'{name}.fa')
                all_data.append(process_file(input_file_path, name))
    df = pd.concat(all_data, ignore_index=True)
    df = df[df['Mutation Info'] != ""]
    df = df.sort_values('Overall Confidence', ascending=False).drop_duplicates('Mutation Info').sort_index()
    # df.to_excel(output_file, index=False)
    df.to_csv(output_file, index=False)

def runbashcmd(cmd, test=False, logf=None):
    print('Executing command:', cmd)
    if test:
        print(cmd)
    err = subprocess.call(cmd, shell=True, stdout=logf, stderr=subprocess.STDOUT)
    if err != 0:
        print('Bash command error: {}\n{}\n'.format(err, cmd))
        sys.exit(1)

def excute_LigandMPNN(input_path, out_dir_path, params_dict, value, env, excute_file):
    cmd = "{} {} --model_type {} --seed {} --pdb_path {} --out_folder {} --ligand_mpnn_cutoff_for_score {} --batch_size {} --number_of_batches {} --redesigned_residues {}".format(
        env + '/bin/python', excute_file, params_dict['model_type'], params_dict['seed'],
        input_path, out_dir_path, params_dict['ligand_mpnn_cutoff_for_score'],
        params_dict['batch_size'], params_dict['number_of_batches'], value
    )
    runbashcmd(cmd)
    return 0

def execute_process(args):
    process_id, input_path, out_dir_path, params_dict, value,  env, excute_file = args
    print(f"Process {process_id} started")
    result = excute_LigandMPNN(input_path, out_dir_path, params_dict, value, env, excute_file)
    print(f"Process {process_id} finished")
    return process_id, result

def batch_run(residues, input_path, out_dir, params_dict, max_workers, env, excute_file):
    pool = multiprocessing.Pool(processes=max_workers)
    results = []
    design_result_dir = os.path.join(out_dir, "design_result")
    if not os.path.exists(design_result_dir):
        os.makedirs(design_result_dir)
    tasks = []
    for value in residues:
        out_dir_path = os.path.join(design_result_dir, f"redesign{value}")
        tasks.append((value, input_path, out_dir_path, params_dict, value, env, excute_file))
    for result in pool.imap_unordered(execute_process, tasks):
        results.append(result)
    pool.close()
    pool.join()
    gc.collect()
    print("All processes have finished")

def zip_directory(directory_path, output_zip):
    if not os.path.exists(directory_path):
        raise ValueError(f"The directory {directory_path} does not exist")
    shutil.make_archive(output_zip, 'zip', directory_path)
    print(f"Directory '{directory_path}' has been zipped into '{output_zip}.zip'")

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
    

def read_residues_from_pdb( site_pdb_path, out_dir ):   

    os.makedirs(out_dir, exist_ok=True)
    site_csv_file = pdb_to_csv( site_pdb_path, out_dir )

    return site_csv_file


def run(site_csv_file, input_path, out_dir, env, excute_file):

    os.makedirs(out_dir, exist_ok=True)
    max_workers = 4
    s0 = time.time()
    task1_df = pd.read_csv(site_csv_file)
    task1_df['pos'] = task1_df['pos'].astype(str)
    task1_df["chain_pos"] = task1_df['chain'] + task1_df['pos']
    residues = list(set(task1_df["chain_pos"]))

    # if len(residues) > 100:
    #     s1 = time.time()
    #     print("compute time cost: %f s"%(s1-s0))
    #     write_num(out_dir,s1-s0,len(residues))
    # else :
        #
    sublists = list(zip(*[iter(residues)]*2))

    name = os.path.splitext(os.path.basename(input_path))[0]
       
    params_dict = {
            "seed": 111,
            "model_type": "ligand_mpnn",
            "ligand_mpnn_cutoff_for_score": "8.0",
            "batch_size": 1,
            "number_of_batches": 20
        }
    for sub in sublists: 
        batch_run(sub, input_path, out_dir, params_dict, max_workers, env, excute_file)  
    final_result_dir = os.path.join(out_dir, "final_result") 
    if not os.path.exists(final_result_dir):
        os.makedirs(final_result_dir)

    process_folder(out_dir, os.path.join(final_result_dir, 'LigandMPNN-overall-output.csv'), name)
    # runbashcmd(f"cp {site_pdb_path} {out_dir}", test=False)
    runbashcmd(f"cp {site_csv_file} {out_dir}", test=False)
    runbashcmd(f"cp {input_path} {out_dir}", test=False)
    zip_directory(out_dir, output_zip=os.path.join(out_dir, "result")) 

    s1 = time.time()
    write_num(out_dir,s1-s0,len(residues))
    print("compute time cost: %f s"%(s1-s0))
        

def main():
    
    parser = argparse.ArgumentParser(description='ProteinMPNN')
    parser.add_argument('--site_pdb', default="", help='site pdb file')
    parser.add_argument('--site_csv', default="", help='site csv file')

    parser.add_argument('--complex_pdb', default="./inputs/BsYjiC.pdb", help='complex pdb file')
    parser.add_argument('--out_dir', default='./outputs/BsYjiC/', help='output dir')
    parser.add_argument('--env', default='/hpcfs/fhome/yangchh/software/anaconda3/envs/ligandmpnn_env', help='env')
    parser.add_argument('--ligandpnn', default='LigandMPNN', help='LigandMPNN path')
    
    args = parser.parse_args()
    if args.site_pdb!="":
        site_csv_file = read_residues_from_pdb( args.site_pdb,  args.out_dir)
        run(site_csv_file, args.complex_pdb, args.out_dir, args.env, args.ligandpnn)
    elif args.site_csv!="":
        run(args.site_csv, args.complex_pdb, args.out_dir, args.env, args.ligandpnn)

if __name__ == "__main__":
    # input1
    # python run_super_hpc.py --site_pdb /hpcfs/fhome/yangchh/ai/ProteinMPNN/1-4CL-1/1-8A1.pdb --complex_pdb /hpcfs/fhome/yangchh/ai/ProteinMPNN/1-4CL-1/1-4CL-1.pdb --out_dir /hpcfs/fhome/yangchh/ai/ProteinMPNN/1-4CL-3 --ligandpnn /hpcfs/fhome/yangchh/ai/ProteinMPNN/LigandMPNN/run.py
    # input2
    # python run_super_hpc.py --site_csv /hpcfs/fhome/yangchh/ai/ProteinMPNN/1-4CL-2/1-8A2.csv   --complex_pdb /hpcfs/fhome/yangchh/ai/ProteinMPNN/1-4CL-1/1-4CL-1.pdb --out_dir /hpcfs/fhome/yangchh/ai/ProteinMPNN/1-4CL-3 --ligandpnn /hpcfs/fhome/yangchh/ai/ProteinMPNN/LigandMPNN/run.py
    main()
