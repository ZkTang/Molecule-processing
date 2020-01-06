import pandas as pd

#load data
df_protein_train = pd.read_csv('F:/ComData/df_protein_train.csv')
df_protein_test = pd.read_csv('F:/ComData/df_protein_test.csv')
protein_merge = pd.concat([df_protein_train, df_protein_test])

#Extincton coefficient
def Ext(protein):
    amino_acid_avail = 'YWC'
    ext_y = 1490
    ext_w = 5500
    ext_c = 125 * 2
    count_aa = {}
    for i in amino_acid_avail:
        count_aa[i] = protein.count(i)
    extinction_coefficient = count_aa['Y'] * ext_y + count_aa['W'] * ext_w + count_aa['C'] * ext_c
    return extinction_coefficient

#output
result = pd.DataFrame([Ext(i.upper()) for i in protein_merge['Sequence']])
result.columns = ['aliphatic_index']
result['aliphatic_index'] = result.aliphatic_index.astype(float)
result['Protein_ID'] = pd.Series(protein_merge['Protein_ID'].values)
result.to_csv('F:/output/Ext_result.csv', sep=',', header=True, index=False)