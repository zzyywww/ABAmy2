#!/usr/bin/env python
#_*_coding:utf-8_*_

from antiberty import AntiBERTyRunner
import torch
import numpy as np
from Bio import SeqIO
import pandas as pd
import sys
import pickle
import math
import time

print('\n' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())+'\n')
input = sys.argv[1]#sequence fasta file
output = sys.argv[2]#resultfile
 

def fasta2list(fasta): 
    name = []
    sequences = []
    for seq_record in SeqIO.parse(fasta, "fasta"):
        name.append(''.join(seq_record.id))
        sequences.append(''.join(seq_record.seq))
    print(f"############# Read {fasta} with {len(sequences)} sequences #############\n")
    return name,sequences

def batch_loader(data, batch_size):
    num_samples = len(data)
    for i in range(0, num_samples, batch_size):
        end_idx = min(i + batch_size, num_samples)
        yield i, end_idx, data[i:end_idx]

def antiberty_emb(sequence_list):#res_code=False，待补充
    antiberty = AntiBERTyRunner()
    batch_size = 100
    n_seqs = len(sequence_list)
    dim = 512
    embeddings = torch.empty((n_seqs, dim))
    
    #当序列的数量小于batch_size时，全部一起处理
    if n_seqs <= batch_size:
        x = antiberty.embed(sequence_list)
        x_mean = [a.mean(axis = 0) for a in x]
        embeddings = torch.stack(x_mean)
    # Use the batch_loader function to iterate through batches
    else:
        n_batches = math.ceil(n_seqs / batch_size)
        i = 1
        for start, end, batch in batch_loader(sequence_list, batch_size):
            # print(start,end)
            # print(f'Batch {i}/{n_batches}\n')
            x = antiberty.embed(batch)
            x_mean = [a.mean(axis = 0) for a in x]
            embeddings[start:end] = torch.stack(x_mean)
            i += 1
    return embeddings


with open( './model/SVM.pkl','rb') as f:
    model = pickle.load(f)

name_list,sequence_list = fasta2list(input)
seq_len = [len(i) for i in sequence_list]

print("############# Computing antiberty features  #############\n")
emb = antiberty_emb(sequence_list)

pred_label = model.predict(emb) 
pred_prob = model.predict_proba(emb)[:,1] 


df = pd.DataFrame()
df['No'] = range(1,len(name_list)+1)
df['Name'] = name_list
df['length'] = seq_len
df['Probability'] = pred_prob
df['Result'] = pred_label

df.to_csv(output,index=False,sep='\t')

print(f"############# The predict results are saved in {output} #############\n")
print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))