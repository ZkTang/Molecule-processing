
import pandas as pd

#load data
df_protein_train = pd.read_csv('F:/ComData/df_protein_train.csv')
df_protein_test = pd.read_csv('F:/ComData/df_protein_test.csv')
protein_merge = pd.concat([df_protein_train, df_protein_test])

# GRAVY coefficient
def GRAVY(protein):
    amino_acid_avail = "ACDEFGHIKLMNPQRSTVWY"
    gravy_para = {'A': 1.80, 'C': 2.50, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2,'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3}
    count_aa = {}
    for i in amino_acid_avail:
        count_aa[i] = protein.count(i)
    length = len(protein)
    value_cal = 0
    for key, value in count_aa.items():
        value_cal = gravy_para[key] * value + value_cal
    gravy_result = value_cal / length
    return round(gravy_result,5)

#output
result = pd.DataFrame([GRAVY(i.upper()) for i in protein_merge['Sequence']])
result.columns = ['GRAVY']
result['GRAVY'] = result.GRAVY.astype(float)
result['Protein_ID'] = pd.Series(protein_merge['Protein_ID'].values)
result.to_csv('F:/output/GRAVY_result.csv', sep=',', header=True, index=False)
