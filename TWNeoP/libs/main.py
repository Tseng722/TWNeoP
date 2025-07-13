
import pandas as pd
from cal_scores import *
from os import path
import logging
from logger_config import TWNeoP_Logger
import argparse
import pathlib

# python /work1791/cindy2270/web/TSA_immuno_prediction/tsa_immuno_prediction/libs/main.py -i /work1791/cindy2270/web/TSA_immuno_prediction/test/myfilecopy3.txt -o /work1791/cindy2270/web/TSA_immuno_prediction/test/test_file/ -n 02
if __name__ == "__main__":
    usage = "%(prog)s <-i --input_file> <-o --output_dir> <-n --name>"
    parser = argparse.ArgumentParser(usage=usage, description='TSA_immunogencity_prediction')
    parser.add_argument('-i', '--input_file', required=True, action='store', help='Provide an input file path. This is required.')
    parser.add_argument('-o', '--output_dir', required=True, action='store', help='Provide an output directory. This is required.')
    parser.add_argument('-n', '--name', required=True, action='store', help='Provide a name. This is required.')
    
    pargs = parser.parse_args()
    FILE_PATH = pargs.input_file
    NAME = pargs.name
    OUTPUT_DIR = pargs.output_dir + NAME
# check file
    if not is_path_exist(OUTPUT_DIR): makedirs(OUTPUT_DIR)
    if not is_path_exist(OUTPUT_DIR +'/tmp/'): makedirs(OUTPUT_DIR +'/tmp/')
    output_dir = path.abspath(path.expanduser(OUTPUT_DIR))

    logger = TWNeoP_Logger(output_dir).get_logger()

    logger.log_info('System initialized')

    delimiter = detect_delimiter(FILE_PATH)

    original_df = pd.read_csv(FILE_PATH, sep=delimiter, header=None, names=['Peptide', 'HLA_Type'],engine = 'python')

    logger.log_info('File readed')
    org_len = len(original_df)

    original_df = original_df.drop_duplicates(subset=['Peptide', 'HLA_Type'])
    rmv_len = len(original_df)

    logger.log_info(f'Duplicates removed - {org_len-rmv_len} pMHC have been removed')

    original_df = original_df.sort_values(by=['HLA_Type'], ascending=True)
    logger.log_info('Data sorted')

    original_df['Length'] = original_df['Peptide'].str.len()



    logger.log_info('Start runing pvacbind')
    try:
        pvac_df = pvac(original_df,output_dir)

    except FileNotFoundError:
        logger.log_error(f'See Pvacbind log which under {OUTPUT_DIR}/pvacbind')
    except Exception as e:
        logger.log_error(f'Pvacbind fail with {e}')
    
    pvac_df = pvac_df[['Epitope Seq','HLA Allele','Median IC50 Score','Median Percentile','cterm_7mer_gravy_score','max_7mer_gravy_score','Best Cleavage Position','Best Cleavage Score','Predicted Stability','Half Life','Stability Rank']]


    logger.log_info('Start runing hydro score')
    try:
        hydro_df,prd_hydro = hydro(original_df) # hydro score
    except Exception as e:
        logger.log_error(f'Hydro fail with {e}')



    logger.log_info('Start runing similarity score')
    try:
        similarity_df = similarity(original_df,NAME,output_dir) # similarity score
    except Exception as e:
        logger.log_error(f'Similarity fail with {e}')


    logger.log_info('Start runing bigmhc score')
    try:
        bigmhc_df = bigmhc(original_df,output_dir)
    except Exception as e:
        logger.log_error(f'Bigmhc fail with {e}')



    logger.log_info('Start runing deepHLApan score')
    try:
        deephlapan_df = deepHLApan(original_df,NAME,output_dir)
    except Exception as e:
        logger.log_error(f'deepHLApan fail with {e}')



    final_df = hydro_df.merge(pvac_df,how='right',left_on=['Peptide','HLA_Type'],right_on = ['Epitope Seq','HLA Allele'], indicator=True)
    final_df.drop(columns=['Epitope Seq','HLA Allele','_merge'], inplace=True)
    final_df = pd.merge(final_df,similarity_df,how='outer',right_on='Peptide',left_on='Peptide')
    final_df.drop(columns=['Counts'], inplace=True)
    final_df = pd.merge(final_df,bigmhc_df,how='outer',left_on=['Peptide','HLA_Type'],right_on = ['pep','mhc'], indicator=True)
    final_df.drop(columns=['_merge','pep','mhc'], inplace=True)
    final_df = pd.merge(final_df,deephlapan_df,how='outer',left_on=['Peptide','HLA_Type'],right_on = ['Peptide','Annotation'])
    final_df = pd.merge(final_df,prd_hydro,how='outer',on=['Peptide','HLA_Type'])
    final_df.rename(columns={'Median IC50 Score': 'IC50', 'Median Percentile': 'Percentile'}, inplace=True)


    logger.log_info('Start runing prediction')
    pred_pro = predictor(final_df)
    final_df['Porioritize Score'] = pred_pro

    final_df = iedb_api(final_df)

    final_df = final_df.rename(columns={'HLA_Type': 'HLA Type', 'IEDB_anno': 'Foreignness Anno','hydro_score':'Hydrophobicity','dissimilarity':'Dissimilarity','foreignness_score':'Foreignness Score','BigMHC_IM':'BigMHC IM','binding score':'Binding Score from DeepHLApan','immunogenic score':'DeepHLApan IM'})
    final_df = final_df[['Peptide','HLA Type','Length','IEDB Qualitative','Foreignness Anno','IC50','Percentile','Binding Score from DeepHLApan','Predicted Stability','Half Life','Stability Rank','Best Cleavage Position','Best Cleavage Score','Hydrophobicity','Dissimilarity','Foreignness Score','BigMHC IM','DeepHLApan IM','Porioritize Score']]
            

    final_df = final_df.round(3)   
    final_df.to_csv(output_dir + '/final.csv',index=False)

    logger.log_info('Successful')