import numpy as np
import math
import pandas as pd
from scipy.spatial.distance import cosine


f = open("../data/species_abundance_summary.txt",'r')
disease_list=[]
microbe_list=[]
microbe_in_health={}
microbe_in_diseases={}
line=0
for i in f:
	if line==0:
		line+=1
		pass
	else:
		i=i.strip('\n').split('\t')
		if i[0] not in disease_list and i[0]!='D006262':
			disease_list.append(i[0])
		if i[2] not in microbe_list:
			microbe_list.append(i[2])
		if i[0]=='D006262':
			microbe_in_health[i[2]]=eval(i[5])
		if i[0]!='D006262':
			if i[2] not in microbe_in_diseases:
				microbe_in_diseases[i[2]]={}
				microbe_in_diseases[i[2]][i[0]]=eval(i[5])
			else:
				microbe_in_diseases[i[2]][i[0]]=eval(i[5])



x=np.zeros((len(microbe_list),len(disease_list)))

for microbe in microbe_in_diseases:
	for disease in microbe_in_diseases[microbe]:
		try:
			value=math.log(microbe_in_diseases[microbe][disease]/microbe_in_health[microbe])
			a=float(value)
			row=microbe_list.index(microbe)
			col=disease_list.index(disease)
			x[row][col]=a
		except:
			row=microbe_list.index(microbe)
			col=disease_list.index(disease)
			x[row][col]=0

print(x)

f2 = open("../data/MASI_v1.0_download_microbesInfo.txt",'r',encoding='utf-8')
f3 = open("../data/microbe_in_microphenoDB_id.txt",'r')

select_microbe=[]
for i in f2:
	i=i.strip('\n').split('\t')
	if i[1] in microbe_list and i[1] not in select_microbe:
		select_microbe.append(i[1])
for i in f3:
	i=i.strip('\n').split('\t')
	if i[1] in microbe_list and i[1] not in select_microbe:
		select_microbe.append(i[1])
print(len(select_microbe))

microbe_network=np.zeros((len(select_microbe),len(select_microbe)))
microbe_feature=np.zeros((len(select_microbe),len(select_microbe)))
for i in select_microbe:
	for j in select_microbe:
		# cos_sim = 1-cosine(x[microbe_list.index(i)],x[microbe_list.index(j)])
		s = np.linalg.norm(x[microbe_list.index(i)]) * np.linalg.norm(x[microbe_list.index(j)])
		if s==0:
			cos=1
		else:
			cos=cosine(x[microbe_list.index(i)],x[microbe_list.index(j)])
		if cos < 0:
			cos=cos*-1
		if cos>1:
			cos=1
		if cos<0.5:
			microbe_network[select_microbe.index(i)][select_microbe.index(j)]=1
		microbe_feature[select_microbe.index(i)][select_microbe.index(j)]=1-cos



np.savetxt('A_microbe.txt',microbe_network,fmt="%.6e")
np.savetxt('microbe_similarity.txt',microbe_feature,fmt="%.6e")
f_out=open('microbe_list_in_A.txt','w')
for i in select_microbe:
    f_out.write(i+'\n')
f_out.close()

