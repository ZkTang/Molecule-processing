import pandas as pd

#load data
df_protein_train = pd.read_csv('F:/ComData/df_protein_train.csv')
df_protein_test = pd.read_csv('F:/ComData/df_protein_test.csv')
protein_merge = pd.concat([df_protein_train, df_protein_test])


def Protein_size_pi(protein):
    #protein size
    amino_acid_avail = "ACDEFGHIKLMNPQRSTVWY"
    aa_residue = {'A': 71.0788, 'C': 103.1388, 'D': 115.0886, 'E': 129.1155, 'F': 147.1766, 'G': 57.0519, 'H': 137.1411,'I': 113.1594, 'K': 128.1741, 'L': 113.1594, 'M': 131.1926, 'N': 114.1038, 'P': 97.1167,'Q': 128.1307, 'R': 156.1875, 'S': 87.0782, 'T': 101.1051, 'V': 99.1326, 'W': 186.2132, 'Y': 163.176}
    count_aa ={}
    size = 0
    for i in amino_acid_avail:
        count_aa[i] = protein.count(i)
    for j in amino_acid_avail:
        size = size + count_aa[j] * aa_residue[j]
    size = (size + 18.01524) / 1000
   #Isoelectric point
    def Isoelectric_point(x):
        COOH = "CDEY" #residue with negative charge
        NH2 = "HKR"   #residue with positive charge
        pi = {'C': 9.0, 'D': 4.0, 'E': 4.5, 'H': 6.4, 'K': 10.4, 'R': 12.0, 'Y': 10.0}
        neg = 0
        pos = 0
        for m in COOH:
            neg += (count_aa[m] * (10 ** x)) / (10 ** x + 10 ** pi[m])
        for n in NH2:
            pos += (count_aa[n] * 10 ** count_aa[n]) / (10 ** x + 10 ** pi[n])
        return neg + 10 ** x / (10 ** x + 10 ** 3.2) - 10 ** 8.2 / (10 ** x + 10 ** 8.2) - pos
    #Bisection method
    min_ele = 3.2
    max_ele = 12
    pi_temp = (min_ele + max_ele) / 2
    for i in range(10):
        if Isoelectric_point(pi_temp) > 0:
            max_ele = pi_temp
            pi_temp = (min_ele + max_ele) / 2
        elif Isoelectric_point(pi_temp) < 0:
            min_ele = pi_temp
            x = (min_ele + max_ele) / 2
    return round(size,5), round(pi_temp,5)

#output
result = pd.DataFrame([Protein_size_pi(i.upper()) for i in protein_merge['Sequence']])
result.columns = ['Protein_size','Isoelectric_point']
result['Protein_size'] = result.Protein_size.astype(float)
result['Isoelectric_point'] = result.Isoelectric_point.astype(float)
result['Protein_ID'] = pd.Series(protein_merge['Protein_ID'].values)
result.to_csv('F:/output/Protein_size_pi.csv', sep=',', header=True, index=False)
