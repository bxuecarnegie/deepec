import gzip

from Bio import SeqIO
from keras.models import load_model
import numpy as np
import sys
import copy


MIN_NUM_OF_SEQS = 10
CLEAN_EC = True
MAX_SEQ_LEN = 1000
# Ignore minim sequence length, original length is 10
MIN_SEQ_LEN = 0

AA_ENCODING = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X']
AA_SCORING = \
    ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', '_',
     'U', 'O', 'B', 'Z', 'J', 'X', '*']
AA_LEN = len(AA_ENCODING)


def preprocessing(fasta_file, temp_file):
    max_len = MAX_SEQ_LEN

    with open(fasta_file, 'r') as input_handle, open(temp_file, 'w') as fp:
        for seq_record in SeqIO.parse(input_handle, "fasta"):
            seq_id = seq_record.id
            seq = seq_record.seq
            if len(seq) <= MAX_SEQ_LEN:
                fp.write('>%s\n' % seq_id)
                fp.write('%s\n' % (seq.strip()))
            else:
                for i in range(0, len(seq) - max_len + 1, 100):
                    new_seq_id = '%s_SEPARATED_SEQUENCE_(%s_%s)' % (seq_id, i + 1, i + max_len + 1)
                    new_seq = seq[i:i + max_len]
                    fp.write('>%s\n' % new_seq_id)
                    fp.write('%s\n' % new_seq)
    return


def fill_aa(seq):
    fill_aa_cnt = MAX_SEQ_LEN - len(seq)
    add_aa_seq = '_' * fill_aa_cnt
    new_seq = seq + add_aa_seq
    return new_seq


def score_info():
    aa_list = AA_SCORING
    aa_score_info = {}
    for aa in aa_list:
        for aa2 in aa_list:
            if aa == aa2:
                aa_score_info[(aa, aa2)] = 1.0
                aa_score_info[(aa2, aa)] = 1.0
            else:
                aa_score_info[(aa, aa2)] = 0.0
                aa_score_info[(aa2, aa)] = 0.0
    return aa_score_info


def one_hot_encoding(seq, aa_score_info):
    data = np.zeros((MAX_SEQ_LEN, AA_LEN), dtype=np.float32)
    aa_list = AA_ENCODING
    for i, aa in enumerate(seq):
        for j, aa2 in enumerate(aa_list):
            data[i, j] = aa_score_info[aa, aa2]
    return data


def run_one_hot_encoding(fasta_file, temp_file):
    aa_score_info = score_info()
    feature_names = ['Feature%s' % i for i in range(1, MAX_SEQ_LEN * AA_LEN + 1)]

    # Using gzip due to the size of output
    with open(fasta_file, 'r') as input_handle, gzip.open(temp_file, 'wt') as fp:
        fp.write('%s\n' % (','.join(['ID'] + feature_names)))

        for seq_record in SeqIO.parse(input_handle, "fasta"):
            try:
                seq_id = seq_record.id
                seq = seq_record.seq

                if MIN_SEQ_LEN <= len(seq) <= MAX_SEQ_LEN:
                    if len(seq) < MAX_SEQ_LEN:
                        seq = fill_aa(seq)
                    encoded_vector = one_hot_encoding(seq, aa_score_info)
                    flatten_encoded_vector = encoded_vector.flatten()
                    flatten_encoded_vector_str = [str(each_val) for each_val in flatten_encoded_vector]
                    fp.write('%s\n' % (','.join([seq_id] + flatten_encoded_vector_str)))
            except (TypeError, ValueError):
                pass


# Merging predicts
def predict_dl(df, output_file, deep_ec_model, multi_label_binarizer=None, threshold=0.5):
    if sys.version_info > (3, 0):
        import _pickle as cPickle
    else:
        import cPickle

    seq_ids = list(df.index)
    X_temp = df.values
    new_X = []
    for i in range(len(X_temp)):
        temp = np.reshape(X_temp[i], (MAX_SEQ_LEN, AA_LEN))
        new_X.append(temp)

    X = np.asarray(new_X)
    X = X.reshape(X.shape[0], MAX_SEQ_LEN, AA_LEN, 1)

    model = load_model(str(deep_ec_model))

    y_predicted = model.predict(X)
    enzyme_list = []

    if multi_label_binarizer is None:
        with open(output_file, 'w') as fp:
            fp.write('Query ID\tPredicted class\tDNN activity\n')
            for i in range(len(y_predicted)):
                score = y_predicted[i][1]
                if y_predicted[i][1] > 0.5:
                    enzyme_list.append(seq_ids[i])
                    fp.write('%s\t%s\t%s\n' % (seq_ids[i], 'Enzyme', score))
                else:
                    fp.write('%s\t%s\t%s\n' % (seq_ids[i], 'Non-enzyme', 1 - score))
    else:
        original_y_predicted = copy.deepcopy(y_predicted)

        y_predicted[y_predicted >= threshold] = 1
        y_predicted[y_predicted < threshold] = 0

        with open(multi_label_binarizer, 'rb') as fid:
            lb = cPickle.load(fid)
        # y_predicted_results = lb.inverse_transform(y_predicted)

        with open(output_file, 'w') as fp:
            fp.write('Query ID\tPredicted EC number\tDNN activity\n')
            for i in range(len(y_predicted)):
                each_y_predicted = copy.deepcopy(y_predicted[i])
                predicted_ddi_score = original_y_predicted[i]
                target_index = np.flatnonzero(each_y_predicted == 1)
                if len(target_index) > 0:
                    for each_idx in target_index:
                        each_y_predicted = copy.deepcopy(y_predicted[i])
                        each_y_predicted[:] = 0
                        each_y_predicted[each_idx] = 1

                        each_y_predicted = np.array([list(each_y_predicted)])
                        y_transformed = lb.inverse_transform(each_y_predicted)

                        score = predicted_ddi_score[each_idx]
                        ec_number = y_transformed[0][0]
                        enzyme_list.append(seq_ids[i])
                        fp.write('%s\t%s\t%s\n' % (seq_ids[i], ec_number, score))

                else:
                    fp.write('%s\t%s\t%s\n' % (seq_ids[i], 'EC number not predicted', 'N/A'))
    return enzyme_list
