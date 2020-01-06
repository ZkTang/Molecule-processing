
import pandas as pd

#load data
df_protein_train = pd.read_csv('F:/ComData/df_protein_train.csv')
df_protein_test = pd.read_csv('F:/ComData/df_protein_test.csv')
protein_merge = pd.concat([df_protein_train, df_protein_test])

#aliphatic index
def aliphatic_index(protein):
    amino_acid_avail = 'AVIL'
    count_aa = {}
    for i in amino_acid_avail:
        count_aa[i] = protein.count(i)
    length = len(protein)
    ai = 100 * (count_aa['A'] / length + 2.9 * count_aa['V'] / length + 3.9 * (count_aa['I'] / length + count_aa['L'] / length))
    return round(ai, 5)

#output
result = pd.DataFrame([aliphatic_index(i.upper()) for i in protein_merge['Sequence']])
result.columns = ['aliphatic_index']
result['aliphatic_index'] = result.aliphatic_index.astype(float)
result['Protein_ID'] = pd.Series(protein_merge['Protein_ID'].values)
result.to_csv('F:/output/aliphatic_index_result.csv', sep=',', header=True, index=False)
