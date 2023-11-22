import re
import subprocess

BLASTP_THREADS = 1


def predict_ec(ref_db_file, target_fasta_file, output_file, output_dir):
    # uniprot_ec_info = {}
    # pfam_uniprot_id_info = {}

    blastp_result = output_dir + '/blastp_result_temp.txt'
    run_blastp(target_fasta_file, blastp_result, ref_db_file)
    try:
        seq_based_ec_prediction_results = read_best_blast_result(blastp_result)

        with open(output_file, 'w') as fp:
            for each_query_id in seq_based_ec_prediction_results:
                all_ec_numbers = seq_based_ec_prediction_results[each_query_id][0].split('|')
                for each_ec_number in sorted(all_ec_numbers):
                    fp.write('%s\t%s\n' % (each_query_id, each_ec_number))
    except OSError:
        pass


def run_blastp(target_fasta, blastp_result, db_dir):
    threads = BLASTP_THREADS
    subprocess.call(
        "diamond blastp -d %s -q %s -o %s --threads %s --id 50 --outfmt 6 qseqid sseqid evalue score qlen slen length "
        "pident" % (
            db_dir, target_fasta, blastp_result, threads), shell=True, stderr=subprocess.STDOUT)


def read_best_blast_result(blastp_result):
    query_db_set_info = {}
    with open(blastp_result, 'r') as fp:
        for line in fp:
            sptlist = line.strip().split('\t')
            query = sptlist[0].strip()
            target = sptlist[1].strip()

            query_id = re.split(r'[|\s]', query)[0]
            ec_numbers = set([ec for ec in re.split(r'[|\s]', target)[1:] if ec != ''])
            score = float(sptlist[3].strip())
            qlen = sptlist[4].strip()
            length = sptlist[6].strip()
            length = float(length)

            coverage = length / float(qlen) * 100
            if coverage >= 75:
                if query_id not in query_db_set_info:
                    query_db_set_info[query_id] = [ec_numbers, score]
                else:
                    p_score = query_db_set_info[query_id][1]
                    # Previously doesn't take into account of multiple best results
                    if score > p_score:
                        query_db_set_info[query_id] = [ec_numbers, score]
                    elif score == p_score:
                        cur_ec_numbers = set(query_db_set_info[query_id][0])
                        cur_ec_numbers.update(ec_numbers)
                        query_db_set_info[query_id] = [cur_ec_numbers, score]
    output_dict = {}
    for query_id in query_db_set_info:
        output_dict[query_id] = ['|'.join(sorted(query_db_set_info[query_id][0])), query_db_set_info[query_id][1]]
    return output_dict
