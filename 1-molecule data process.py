import numpy as np
import pandas as pd

# molcule data
df_molecule = pd.read_csv('F:\ComData\df_molecule.csv')

#fingerprint
molecule_fingerprint = (df_molecule['Fingerprint'].str.split(',')).tolist()
molecule_fingerprint = pd.DataFrame(molecule_fingerprint)
molecule_fingerprint.columns = ['finger_{}'.format(i) for i in range(1, 168)]
df_molecule = df_molecule.drop(['Fingerprint'], axis=1)
#print(molecule_fingerprint)

#feature & missing value process
df_molecule_feat = df_molecule.drop(['Molecule_ID'], axis=1)
df_molecule_feat['cyp_3a4'] = df_molecule_feat['cyp_3a4'].fillna(df_molecule_feat['cyp_3a4'].mean())
df_molecule_feat['cyp_2c9'] = df_molecule_feat['cyp_2c9'].fillna(df_molecule_feat['cyp_2c9'].mean())
df_molecule_feat['cyp_2d6'] = df_molecule_feat['cyp_2d6'].fillna(df_molecule_feat['cyp_2d6'].mean())
df_molecule_feat['ames_toxicity'] = df_molecule_feat['ames_toxicity'].fillna(df_molecule_feat['ames_toxicity'].quantile(0.75))
df_molecule_feat['fathead_minnow_toxicity'] = df_molecule_feat['fathead_minnow_toxicity'].fillna(df_molecule_feat['fathead_minnow_toxicity'].mode()[0])
df_molecule_feat['tetrahymena_pyriformis_toxicity'] = df_molecule_feat['tetrahymena_pyriformis_toxicity'].fillna(df_molecule_feat['tetrahymena_pyriformis_toxicity'].mode()[0])
df_molecule_feat['honey_bee'] = df_molecule_feat['honey_bee'].fillna(df_molecule_feat['honey_bee'].quantile(0.75))
df_molecule_feat['logP'] = df_molecule_feat['logP'].fillna(df_molecule_feat['logP'].mean())
df_molecule_feat['CLtotal'] = df_molecule_feat['CLtotal'].fillna(df_molecule_feat['CLtotal'].mean())
df_molecule_feat['hia'] = df_molecule_feat['hia'].fillna(df_molecule_feat['hia'].mode()[0])
df_molecule_feat['biodegradation'] = df_molecule_feat['biodegradation'].fillna(df_molecule_feat['biodegradation'].mode()[0])
df_molecule_feat['Vdd'] = df_molecule_feat['Vdd'].fillna(df_molecule_feat['Vdd'].mode()[0])
df_molecule_feat['p_glycoprotein_inhibition'] = df_molecule_feat['p_glycoprotein_inhibition'].fillna(df_molecule_feat['p_glycoprotein_inhibition'].mean())
df_molecule_feat['NOAEL'] = df_molecule_feat['NOAEL'].fillna(df_molecule_feat['NOAEL'].mean())
df_molecule_feat['bbb'] = df_molecule_feat['bbb'].fillna(df_molecule_feat['bbb'].mode()[0])

#output
df_molecule_feat['Molecule_ID'] = df_molecule['Molecule_ID']
df_molecule_all = pd.concat([df_molecule_feat, molecule_fingerprint],axis=1)
df_molecule_all.to_csv('F:/output/molecule.csv', index=False)










