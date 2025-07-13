from os import path, makedirs, getcwd
import glob
import os
import pandas as pd
import subprocess
import math
import logging
from joblib import load
import pathlib


PACKAGE_DIR = pathlib.Path(__file__).resolve().parents[1] 

def is_path_exist(dir, error_msg=False):
    if dir == None: return False
    elif path.exists(dir): return True
    else: return False

def detect_delimiter(file_path):
    try:
        df_comma = pd.read_csv(io.StringIO(file_path), sep=',', header=None, nrows=1)
        if len(df_comma.columns) > 1:
            return ','
    except Exception:
        pass
    try:
        df_tab = pd.read_csv(io.StringIO(file_path), sep='\t', header=None, nrows=1)
        if len(df_tab.columns) > 1:
            return '\t'
    except Exception:
        pass
    try:
        df_comma = pd.read_csv(io.StringIO(file_path), sep=' ', header=None, nrows=1)
        if len(df_comma.columns) > 1:
            return ' '
    except Exception:
        pass

def run_my_subprocess(command):
    try:
        subprocess.run(command, shell=True, check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
      
    except subprocess.CalledProcessError as e:
        print(f"command: {command}")
        print(f"Error: {e}")


def pvac(df,output_dir):
    output_dir = os.path.join(output_dir, 'pvacbind')
    if not is_path_exist(output_dir): makedirs(output_dir)
    
    for i in range(len(df)):
        hla = df.at[i,'HLA_Type']
        pep = df.at[i,'Peptide']
        l = len(pep)
        file = open(f'{output_dir}/{hla}_{l}.fasta','a')
        file.write(">" + str(pep) +'_' + hla +'_' + '\n'  + str(pep)+ '\n')
    file.close()
    pattern = os.path.join(output_dir, '*.fasta')
    fasta_files = glob.glob(pattern)
    final_df = pd.DataFrame()
    for file in fasta_files:
        file_name = os.path.splitext(os.path.basename(file))[0]
        hla = file_name.split('_')[0]
        mer = file_name.split('_')[1]
        try:
            cmd = f'pvacbind run  {file}  {file_name}  {hla}  NetMHCpan NetMHC  {output_dir}/{file_name}/  -e1 {mer}   -k  -t 48 -b 100000  --aggregate-inclusion-binding-threshold 100000  --net-chop-method cterm  --netmhc-stab   --iedb-install-directory {PACKAGE_DIR}/supporting_data/iedb '   
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output = True, 
                text = True 
            )
            o = open(f'{output_dir}/{file_name}/out.txt', 'w')
            e = open(f'{output_dir}/{file_name}/error.txt', 'w')
            o.write(result.stdout)
            e.write(result.stderr)
            o.close()
            e.close()
        
        except Exception as e:
            print(f"An exception occurred: {str(e)}")

        tmp_dir = f'{output_dir}/{file_name}/MHC_Class_I/{file_name}.filtered.tsv'
        tmp_df = pd.read_csv(tmp_dir,sep="\t")
        final_df = pd.concat([final_df,tmp_df],ignore_index=True)

    return final_df


def hydro_vector(pep):
	hydro_score={"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-0.4,"H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,"P":-1.6,"Q":-3.5,"R":-4.5,"S":-0.8,"T":-0.7,"V":4.2,"W":-0.9,"Y":-1.3}
	hydrophobicity_vector=[]
	pep_list=list(pep)
	pep_len=len(pep_list)
	for pep in pep_list:
		hydrophobicity_vector.append(hydro_score[pep.upper()])
	return hydrophobicity_vector
def calculate_hydro(peptide,hla,mer,df_weight):
	f = (df_weight['HLA allele'] ==hla)
	df_tmp = df_weight.loc[f].reset_index(drop=True)
	if df_tmp.empty:
		return 'NA','NA'
    
	hydro_final =[]
	hydro_list = hydro_vector(peptide)
	hydro_score_total = 0
	try:
		for i in range(mer): 
			hydro_score_single = (-math.log(df_tmp.iat[0,i+1],10)*hydro_list[i])
			hydro_final.append(hydro_score_single)
			hydro_score_total += hydro_score_single
		hydro_score_avg = hydro_score_total/mer
	except:
		hydro_score_avg = 'NA','NA'
	return hydro_score_avg,hydro_final

def hydro(df):
    df_9_mer = pd.read_excel("/CMU_TSA/cindy2270/IEDB/score_TCR/abg2200_Data_file_S2.xlsx",sheet_name=0,index_col=None)
    df_8_mer = pd.read_excel("/CMU_TSA/cindy2270/IEDB/score_TCR/abg2200_Data_file_S2.xlsx",sheet_name=1,index_col=None)
    df_10_mer = pd.read_excel("/CMU_TSA/cindy2270/IEDB/score_TCR/abg2200_Data_file_S2.xlsx",sheet_name=2,index_col=None)
    df_11_mer = pd.read_excel("/CMU_TSA/cindy2270/IEDB/score_TCR/abg2200_Data_file_S2.xlsx",sheet_name=3,index_col=None)
    df_hydro = pd.DataFrame()
    for i in range(len(df)):
        pep = df.at[i,'Peptide']
        hla = df.at[i,'HLA_Type']
        mer = df.at[i,'Length']
        df_hydro.at[i,'Peptide'] = pep 
        df_hydro.at[i,'HLA_Type'] = hla
        if mer == 8:
            df_weight = df_8_mer
        elif mer == 9:
            df_weight = df_9_mer
        elif mer == 10:
            df_weight = df_10_mer
        elif mer == 11:
            df_weight = df_11_mer
        else:
            df.at[i,'hydro_score'] = 0
            df_hydro.at[i,'hydro_avg_score'] = 0
            continue
        hydro_score, hydro_final= calculate_hydro(pep,hla,mer,df_weight)
        df.at[i,'hydro_score'] = hydro_score
        df_hydro.at[i,'hydro_avg_score'] = hydro_score
        if hydro_final=='NA':
            continue
        else:
            for j in range(1,mer+1):
                df_hydro.at[i,f'Position {j}'] = hydro_final[j-1]

    return df,df_hydro



def similarity(df,sample_name,output):
    # container_id = "sleepy_dirac"
    container_id = "similarity"
    
    # r_script = f"Rscript /root/local179/test_docker.R"
    r_script = f"Rscript /root/similarity_docker.R"
    output_txt = output +f'/tmp/{sample_name}.txt'
    
    dfg = df.groupby(['Peptide']).size().reset_index(name="Counts")
    dfg['Peptide'].to_csv(output_txt, header=False, index=False)
    
    
    command = [
        [f"docker run -v {output}/tmp:/root/local179/ -it -d --name {container_id} andrewrech/antigen.garnish:2.3.0 /bin/bash"],
        [f"docker cp {PACKAGE_DIR}/supporting_data/similarity/netMHC-4.0a.Linux.tar.gz {container_id}:/netMHC-4.0a.Linux.tar.gz"],
        [f"docker cp {PACKAGE_DIR}/supporting_data/similarity/netMHCII-2.3.Linux.tar.gz {container_id}:/netMHCII-2.3.Linux.tar.gz"],
        [f"docker cp {PACKAGE_DIR}/supporting_data/similarity/netMHCpan-4.1b.Linux.tar.gz {container_id}:netMHCpan-4.1b.Linux.tar.gz"],
        [f"docker cp {PACKAGE_DIR}/supporting_data/similarity/netMHCIIpan-4.0.Linux.tar.gz {container_id}:netMHCIIpan-4.0.Linux.tar.gz"],
        [f"docker exec {container_id} config_netMHC.sh"],
        [f"docker cp {PACKAGE_DIR}/libs/similarity_docker.R {container_id}:/root/similarity_docker.R"],
        # [f"docker cp {output_txt} {container_id}:/root/{sample_name}.txt"],
        [f'docker start {container_id}'],
        [f'docker exec -e PATH="/root/antigen.garnish/netMHC/netMHC-4.0:/root/antigen.garnish/netMHC/netMHCII-2.3:/root/antigen.garnish/netMHC/netMHCIIpan-4.0/:/root/antigen.garnish/netMHC/netMHCpan-4.1:/root/antigen.garnish/ncbi-blast-2.10.1+-src/c++/ReleaseMT/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" -e AG_DATA_DIR="/root/antigen.garnish" {container_id} {r_script} {sample_name}'],
        [f'docker stop {container_id}'],
        [f'docker rm {container_id}']
    ]
    for cmd in command:
        run_my_subprocess(cmd)
    
    df_f = pd.read_table(output +f'/tmp/f_{sample_name}.txt')
    df_d = pd.read_table(output +f'/tmp/d_{sample_name}.txt')
    df_mf = pd.merge(dfg,df_f,how='outer',right_on='nmer',left_on='Peptide')
    df_mf.drop('nmer', axis=1, inplace=True)
    df_md = pd.merge(df_mf,df_d,how='outer',right_on='nmer',left_on='Peptide')
    df_md.drop('nmer', axis=1, inplace=True)


    if 'IEDB_anno' not in df_md.columns:
        df_md['IEDB_anno'] = ''
    
    return df_md


def bigmhc(df,output_dir):
    file_path = output_dir + '/tmp/bigmhc_tmp.csv'
    df.to_csv(file_path,index =False)
    output_file = output_dir+'/tmp/bigmhc_result.csv'
    command = f'python {PACKAGE_DIR}/libs/bigmhc/src/predict.py -i={file_path} -m=im -a=1 -p=0 -c=1 -o={output_file} -d="cpu"'
    try:
        subprocess.run(command, shell=True, check=True,stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
    df_bigmhc = pd.read_csv(output_file)
    df_bigmhc.drop(columns=['tgt','len'], inplace=True)
    return df_bigmhc

def deepHLApan(df_in,name,output):
    # container_id = "deephla"
    container_id = "deephlapan"
    df = df_in
    df['Annotation'] = df['HLA_Type']
    df['HLA'] = df['HLA_Type'].str.replace('*','')
    df['peptide'] = df['Peptide']
    exclude_hla = ['HLA-B30:02', 'HLA-A33:02', 'HLA-A35:01', 'HLA-A51:01', 'HLA-B07:01']
    for hla in exclude_hla:
        df = df.loc[~df['HLA'].str.contains(hla, na=False)]
    output_csv = output+f'/tmp/deepHLApan_{name}.csv'
    df[['Annotation','HLA','peptide']].to_csv(output_csv,index=False)
    
    commands = [
        [f"docker run -v {output}/tmp:/root/local179 -it -d --name {container_id} biopharm/deephlapan:v1.1 bash"],
        [f"docker start {container_id}"],
        [f'docker exec {container_id} deephlapan -F /root/local179/deepHLApan_{name}.csv -O /root/local179'],
        [f'docker stop {container_id}'],
        [f'docker rm {container_id}']
    ]

    for cmd in commands:
        run_my_subprocess(cmd)

    df_dhp = pd.read_csv(f'{output}/tmp/deepHLApan_{name}_predicted_result.csv')
    df_dhp.drop(columns='HLA',inplace=True)
    
    return df_dhp


def predictor(df):
    chk_len = len(df)
    feature_cols = ['Median IC50 Score', 'Median Percentile','Best Cleavage Position','foreignness_score','Best Cleavage Score',
        'Predicted Stability','dissimilarity',
       'Half Life', 'Stability Rank','hydro_avg_score', 'Position 1',
       'Position 2', 'Position 3', 'Position 4', 'Position 5', 'Position 6',
       'Position 7', 'Position 8','Position 9','Position 10','Position 11','binding score', 'BigMHC_IM','immunogenic score']
    df = df.reindex(columns=feature_cols, fill_value=0)
    df.fillna(0, inplace=True)
    
    loaded_model = load(f'{PACKAGE_DIR}/supporting_data/predictor/best-089.joblib')
    pred_pro = loaded_model.predict_proba(df)[:, 1]
    if len(pred_pro) != chk_len:
        return False
    return pred_pro



def iedb_api(df):
    df_iedb = pd.DataFrame()
    for i in range(len(df)):
        df_tmp = pd.DataFrame()
        tumor = df.at[i,'Peptide']
        hla = df.at[i,'HLA_Type']
        try:
            r = requests.get(f'https://query-api.iedb.org/mhc_search?linear_sequence=eq.{tumor}&mhc_restriction=eq.{hla}&select=linear_sequence%2Cstructure_type%2Csource_organism_name%2Celution_id%2Cmhc_restriction%2Cqualitative_measure%2Chost_organism_name%2Cmhc_allele_name')
            dict_obj = json.loads(r.text)
            if not dict_obj:
                df_tmp.at[0,'IEDB Qualitative'] = ''
                df_tmp.at[0,'Peptide'] = tumor
                df_tmp.at[0,'HLA_Type'] = hla
            else:
                pos_cpunt = 0
                neg_vount = 0
                for j,v in enumerate(dict_obj):
                    if  'osi' in v['qualitative_measure']:
                        pos_cpunt += 1
                    elif 'ega' in v['qualitative_measure']:
                        neg_vount += 1
                df_tmp.at[0,'Peptide'] = tumor
                df_tmp.at[0,'HLA_Type'] = hla 
                df_tmp.at[0,'IEDB Qualitative'] = f'Positive : {pos_cpunt} ; Negative : {neg_vount}'
        except:
            df_tmp.at[0,'IEDB Qualitative'] = ''
            df_tmp.at[0,'Peptide'] = tumor
            df_tmp.at[0,'HLA_Type'] = hla
        df_iedb = pd.concat([df_iedb,df_tmp],axis = 0)
    df = pd.merge(df,df_iedb,how='outer',right_on=['Peptide','HLA_Type'],left_on=['Peptide','HLA_Type'])
    return df

def test():
    
    print(package_dir)