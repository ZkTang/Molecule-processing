import pandas as pd
import numpy as np
import re
from gensim.models import Word2Vec

#load data
df_protein_train = pd.read_csv("F:\ComData\df_protein_train.csv")
df_protein_test = pd.read_csv("F:\ComData\df_protein_test.csv")
df_protein = pd.concat([df_protein_train,df_protein_test])
df_protein.Sequence = df_protein.Sequence.apply(lambda x: x.upper())
#exit()

# set parameters & split sequence
embed_size = 128
seq_size = 3
w2v_column = [f'w2v_{i}' for i in range(embed_size)]
sentence_split = []
for item in df_protein.Sequence:
        for t in range(0,seq_size):
            sentence_split.append([word for word in re.findall(r'.{'+str(seq_size)+'}',item[t:])])
#print(sentence_split)

#set model
model = Word2Vec(sentence_split,size = embed_size,window = 4,min_count = 1,negative = 3,sg = 1,sample = 0.001,hs = 1,workers = 4,iter = 15)
w2v_feat = []

#generate feature
i = 0
while i <= len(sentence_split) - seq_size:
        sum_w2v = np.zeros(shape = (embed_size,))
        for j in range(i,i + seq_size):
            for word in sentence_split[j]:
                sum_w2v += model[word]
        w2v_feat.append(sum_w2v)
        i = i + seq_size

#output
w2v_feat = np.vstack(w2v_feat)
df_w2v = pd.DataFrame(w2v_feat,columns=w2v_column)
df_w2v["Protein_ID"] = df_protein["Protein_ID"].values
df_w2v.to_csv('F:/output/w2v_result.csv', sep=',', header=True, index=False)