import numpy as np
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit import DataStructs
import operator

NPs = []
SMIs = []
f1 = open("../data/All_CIDs_matched_SMILES.smi",'r')
for lines in f1:
	lines = lines.strip().split('\t')
	smi = lines[0]
	np_id = lines[1]
	NPs.append(np_id)
	SMIs.append(smi)
f1.close()

f_drug=open('drug_list_in_A.txt','w')
for each in NPs:
	f_drug.write(each+'\n')


S = np.zeros((len(NPs),len(NPs)))

fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(i)) for i in SMIs]

for i in range(len(fps)):
    for j in range(len(fps)):
        S[i][j] = DataStructs.CosineSimilarity(fps[i], fps[j])
print(S)

np.savetxt('drug_similarity.txt',S,fmt="%.6e")
A_drug=np.zeros(shape=(S.shape[0],S.shape[1]))
for row in range(S.shape[0]):
	for col in range(S.shape[1]):
		if S[row][col] > 0.5:
			A_drug[row][col]=1
		else:
			A_drug[row][col]=0
np.savetxt('A_drug.txt',A_drug,fmt="%.6e")

