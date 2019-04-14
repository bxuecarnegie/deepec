from Bio import SeqIO
import gzip
import shutil
import logging
import subprocess

def predict_ec(ref_db_file, target_fasta_file, output_file, output_dir):
    uniprot_ec_info = {}
    pfam_uniprot_id_info = {}

    blastp_result = output_dir+'/blastp_result_temp.txt'
    run_blastp(target_fasta_file, blastp_result, ref_db_file) 
    try:
        seq_based_ec_prediction_results = read_best_blast_result(blastp_result)

        with open(output_file, 'w') as fp:
            for each_query_id in seq_based_ec_prediction_results:
                each_ec_number = seq_based_ec_prediction_results[each_query_id][0]
                fp.write('%s\t%s\n'%(each_query_id, each_ec_number))
    except:
        pass

def run_blastp(target_fasta, blastp_result, db_dir):    
    threads = 1
    subprocess.call("diamond blastp -d %s -q %s -o %s --threads %s --id 50 --outfmt 6 qseqid sseqid evalue score qlen slen length pident"%(db_dir, target_fasta, blastp_result, threads), shell=True, stderr=subprocess.STDOUT)

def read_best_blast_result(blastp_result):
    query_db_set_info = {}
    with open(blastp_result, 'r') as fp:
        for line in fp:
            sptlist = line.strip().split('\t')            
            query_id = sptlist[0].strip()
            db_id = sptlist[1].strip()  
            
            ec_number = db_id.split('|')[1].strip()
            score = float(sptlist[3].strip())
            qlen = sptlist[4].strip()
            length = sptlist[6].strip()
            length = float(length)
            
            coverage = length/float(qlen)*100
            if coverage >= 75:
                if query_id not in query_db_set_info:
                    query_db_set_info[query_id] = [ec_number, score]
                else:
                    p_score = query_db_set_info[query_id][1]
                    if score > p_score:
                        query_db_set_info[query_id] = [ec_number, score]
    return query_db_set_info



