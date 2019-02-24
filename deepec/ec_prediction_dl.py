from Bio import SeqIO
from keras.models import load_model
import numpy as np
import pandas as pd
import sys
import os
import glob
import shutil
import logging
import time
import subprocess
import copy

    
def preprocessing(fasta_file, temp_file):
    l = 1000
    fp = open(temp_file, 'w')
    input_handle = open(fasta_file, "r")
    for seq_record in SeqIO.parse(input_handle, "fasta"):   
        seq_id = seq_record.id
        seq = seq_record.seq
        if len(seq) <= 1000:
            fp.write('>%s\n'%(seq_id))
            fp.write('%s\n'%(seq.strip()))
        else:
            for i in range(0, len(seq)-l+1, 100):
                new_seq_id = '%s_SEPARATED_SEQUENCE_(%s_%s)' % (seq_id, i+1, i+l+1)
                new_seq = seq[i:i+l]
                fp.write('>%s\n'%(new_seq_id))
                fp.write('%s\n'%(new_seq))
    
    input_handle.close()
    fp.close()
    return

def fill_aa(seq):
    fill_aa_cnt = 1000 - len(seq)
    add_aa_seq = '_' * fill_aa_cnt
    new_seq = seq + add_aa_seq
    return new_seq

def score_info():
    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', '_']
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
    data = np.zeros((1000, 21), dtype=np.float32)
    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X']
    for i, aa in enumerate(seq):
        for j, aa2 in enumerate(aa_list):
            data[i, j] = aa_score_info[aa, aa2]
    return data

def run_one_hot_encoding(fasta_file, temp_file):
    aa_score_info = score_info()
    fp = open(temp_file, 'w')
    feature_names = ['Feature%s'%(i) for i in range(1, 21001)] 
    
    fp.write('%s\n'%(','.join(['ID']+feature_names)))
    
    input_handle = open(fasta_file, "r")
    for seq_record in SeqIO.parse(input_handle, "fasta"):   
        try:
            seq_id = seq_record.id
            seq = seq_record.seq
            if len(seq) >= 10 and len(seq) <= 1000:
                if len(seq) < 1000:
                    seq = fill_aa(seq)
                encoded_vector = one_hot_encoding(seq, aa_score_info)
                flatten_encoded_vector = encoded_vector.flatten()
                flatten_encoded_vector_str = [str(each_val) for each_val in flatten_encoded_vector]
                fp.write('%s\n'%(','.join([seq_id]+flatten_encoded_vector_str)))
        except:
            pass
    input_handle.close()
    fp.close()


def predict_target_ec(df, output_file, DeepEC_model):
    seq_ids = list(df.index)
    X_temp = df.values
    new_X = []
    for i in range(len(X_temp)):
        temp = np.reshape(X_temp[i], (1000, 21))
        new_X.append(temp)
    
    X = np.asarray(new_X)
    X = X.reshape(X.shape[0], 1000, 21, 1)
    
    model = load_model(DeepEC_model)
    
    y_predicted = model.predict(X)
    enzyme_list = []
    with open(output_file, 'w') as fp:
        fp.write('Query ID\tPredicted class\tDNN activity\n')
        for i in range(len(y_predicted)): 
            socre = y_predicted[i][1]
            if y_predicted[i][1] > 0.5:
                enzyme_list.append(seq_ids[i])
                fp.write('%s\t%s\t%s\n'%(seq_ids[i], 'Target', socre))
            else:
                fp.write('%s\t%s\t%s\n'%(seq_ids[i], 'Non-target', 1-socre))
    return enzyme_list

def predict_enzyme(df, output_file, DeepEC_model):
    seq_ids = list(df.index)
    X_temp = df.values
    new_X = []
    for i in range(len(X_temp)):
        temp = np.reshape(X_temp[i], (1000, 21))
        new_X.append(temp)
    
    X = np.asarray(new_X)
    X = X.reshape(X.shape[0], 1000, 21, 1)
    
    model = load_model(DeepEC_model)
    
    y_predicted = model.predict(X)
    enzyme_list = []
    with open(output_file, 'w') as fp:
        fp.write('Query ID\tPredicted class\tDNN activity\n')
        for i in range(len(y_predicted)): 
            socre = y_predicted[i][1]
            if y_predicted[i][1] > 0.5:
                enzyme_list.append(seq_ids[i])
                fp.write('%s\t%s\t%s\n'%(seq_ids[i], 'Enzyme', socre))
            else:
                fp.write('%s\t%s\t%s\n'%(seq_ids[i], 'Non-enzyme', 1-socre))
    return enzyme_list

def predict_ec(df, output_file, DeepEC_model, MultiLabelBinarizer, threshold=0.5):
    if (sys.version_info > (3, 0)):
        import _pickle as cPickle
    else:
        import cPickle
    
    seq_ids = list(df.index)
    X_temp = df.values
    new_X = []
    for i in range(len(X_temp)):
        temp = np.reshape(X_temp[i], (1000, 21))
        new_X.append(temp)
    
    X = np.asarray(new_X)
    X = X.reshape(X.shape[0], 1000, 21, 1)
    
    with open(MultiLabelBinarizer, 'rb') as fid:
        lb = cPickle.load(fid)
    
    model = load_model(DeepEC_model)
    
    y_predicted = model.predict(X)
    original_y_predicted = copy.deepcopy(y_predicted)
    
    y_predicted[y_predicted >= threshold] = 1
    y_predicted[y_predicted <threshold] = 0

    y_predicted_results = lb.inverse_transform(y_predicted) 

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
                    fp.write('%s\t%s\t%s\n'%(seq_ids[i], ec_number, score))
            
            else:
                fp.write('%s\t%s\t%s\n'%( seq_ids[i], 'EC number not predicted', 'N/A'))




