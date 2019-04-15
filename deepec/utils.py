import argparse
import copy
import os
import shutil
import pandas as pd
from Bio import SeqIO


def argument_parser(version=None):
    parser = argparse.ArgumentParser()    
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory")
    parser.add_argument('-i', '--fasta_file', required=True, help="Input fasta file") 
    parser.add_argument('-t', '--threshold', default=0.5, help="DNN activity") 
    return parser

def remove_files(output_dir):
    try:
        os.remove(output_dir+'/SwissProt_input.csv')     
    except:
        pass
    try:
        os.remove(output_dir+'/low_seq_candidates.txt')     
    except:
        pass
    try:
        os.remove(output_dir+'/blastp_result_temp.txt')     
    except:
        pass
    try:
        os.remove(output_dir+'/temp_fasta.fa')     
    except:
        pass
    return

def move_files(output_dir):
    try:
        os.mkdir(output_dir+'/log_files/')
    except:
        pass
    
    try:
        shutil.move(output_dir+'/3digit_EC_prediction.txt', output_dir+'/log_files/3digit_EC_prediction.txt')
    except:
        pass
    try:
        shutil.move(output_dir+'/4digit_EC_prediction.txt', output_dir+'/log_files/4digit_EC_prediction.txt')
    except:
        pass
    try:
        shutil.move(output_dir+'/DeepEC_Result_DL.txt', output_dir+'/log_files/DeepEC_Result_DL.txt')
    except:
        pass
    try:
        shutil.move(output_dir+'/Enzyme_prediction.txt', output_dir+'/log_files/Enzyme_prediction.txt')
    except:
        pass
    try:
        shutil.move(output_dir+'/Blastp_result.txt', output_dir+'/log_files/Blastp_result.txt')
    except:
        pass
    return

def read_ec_list(filename):
    ec_list = []
    with open(filename, 'r') as fp:
        for line in fp:
            ec_list.append(line.strip())
    return ec_list

def read_prediction_results(result_file):
    ec_results = {}
    with open(result_file, 'r') as fp:
        fp.readline()
        for line in fp:
            sptlist = line.strip().split('\t')
            query = sptlist[0].strip()
            if '_SEPARATED_SEQUENCE_' in query:
                query = query.split('_SEPARATED_SEQUENCE_')[0].strip()
            predicted_ec = sptlist[1].strip()
            dnn_activity = sptlist[2].strip()
            if predicted_ec != 'EC number not predicted':
                if query not in ec_results:
                    ec_results[query] = [predicted_ec]
                else:
                    ec_results[query].append(predicted_ec)
                    ec_results[query] = list(set(ec_results[query]))

    return ec_results

def get_prediction_results(ec_4digit_results, ec_3digit_results, enzyme_results):
    final_one_digit_ec_results = {}
    for query in ec_4digit_results:
        if query in enzyme_results:
            predicted_ec_list = ec_4digit_results[query]
            for predicted_ec in predicted_ec_list:
                flag3 = False
                if query in ec_3digit_results:
                    for each_digit_ec in ec_3digit_results[query]:
                        if each_digit_ec in predicted_ec:
                            flag3 = True        
                
                if flag3 == True:
                    if query not in final_one_digit_ec_results:
                        final_one_digit_ec_results[query] = [predicted_ec]
                    else:
                        final_one_digit_ec_results[query].append(predicted_ec)
                        final_one_digit_ec_results[query] = list(set(final_one_digit_ec_results[query]))  
                        
    return final_one_digit_ec_results

def merge_dl_prediction_results(output_dir, dl4_output_file, dl3_output_file, enzyme_output_file, fasta_file, target_fasta_file):
    enzyme_results = {}
    with open(enzyme_output_file, 'r') as fp:
        fp.readline()
        for line in fp:
            sptlist = line.strip().split('\t')
            query = sptlist[0].strip()
            if '_SEPARATED_SEQUENCE_' in query:
                query = query.split('_SEPARATED_SEQUENCE_')[0].strip()
                    
            enzyme_info = sptlist[1].strip()
            if 'Non-enzyme' not in enzyme_info:
                enzyme_results[query] = 1
         
    ec_4digit_results = read_prediction_results(dl4_output_file)
    ec_3digit_results = read_prediction_results(dl3_output_file)

    final_one_digit_ec_results = get_prediction_results(ec_4digit_results, ec_3digit_results, enzyme_results)
    
    with open('%s/DeepEC_Result_DL.txt'%(output_dir), 'w') as fp:
        fp.write('Query ID	Predicted EC number\n')
        for each_query in final_one_digit_ec_results:
            for each_predicted_ec in final_one_digit_ec_results[each_query]:
                fp.write('%s\t%s\n'%(each_query, each_predicted_ec))
                
    fp = open(target_fasta_file, 'w')      
    input_handle = open(fasta_file, "r")
    for seq_record in SeqIO.parse(input_handle, "fasta"):   
        seq_id = seq_record.id
        if seq_id in enzyme_results:
            if seq_id not in final_one_digit_ec_results:
                fp.write('>%s\n'%(seq_id))
                fp.write('%s\n'%(seq_record.seq))
    fp.close()
    
def merge_dl_seq_results(output_dir):
    fp2 = open('%s/DeepEC_Result.txt'%(output_dir), 'w')
    fp2.write('Query ID	Predicted EC number\n')
    with open('%s/DeepEC_Result_DL.txt'%(output_dir), 'r') as fp:
        fp.readline()
        for line in fp:
            fp2.write('%s\n'%(line.strip()))
    try:        
        with open('%s/Blastp_result.txt'%(output_dir), 'r') as fp:
            for line in fp:
                fp2.write('%s\n'%(line.strip()))
    except:
        pass
    fp2.close()