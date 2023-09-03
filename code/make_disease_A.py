import re
import numpy as np
from scipy.spatial.distance import cosine
f=open('../data/CTD_diseases_pathways.csv','r')

disease=[]
pathway=[]
dis_pathway={}
for i in f:
	if i[0]!= '#':
		i=i.strip('\n')
		j= re.split(',(?=([^\"]*\"[^\"]*\")*[^\"]*$)',i)
		if j[2] not in disease:
			disease.append(j[2])
		if j[8] not in pathway:
			pathway.append(j[8])
		if j[2] not in dis_pathway:
			dis_pathway[j[2]]=set()
			dis_pathway[j[2]].add(j[8])
		if j[2] in dis_pathway:
			dis_pathway[j[2]].add(j[8])

dis_gene_matrix=np.zeros(shape=(len(disease),len(pathway)))
for each_dis in dis_pathway:
	for each_path in dis_pathway[each_dis]:
		dis_gene_matrix[disease.index(each_dis)][pathway.index(each_path)]=1
# 		print(disease.index(each_dis),pathway.index(each_path))
print(dis_gene_matrix)

dis_cos_matrix=np.zeros(shape=(len(disease),len(disease)))
for row in range(len(disease)):
	for col in range(len(disease)):
		disease1=dis_gene_matrix[row]
		disease2=dis_gene_matrix[col]
		cos_sim = 1-cosine(disease1,disease2)
		dis_cos_matrix[row][col]=cos_sim
print(dis_cos_matrix)

select_disease_id=[]
select_disease=[]
for i in range(len(disease)):
    if sum(dis_gene_matrix[i])!=1:
        select_disease_id.append(i)
        select_disease.append(disease[i])

selected_cos_matrix=np.zeros(shape=(len(select_disease),len(select_disease)))
for row in range(len(select_disease)):
    for col in range(len(select_disease)):
        if dis_cos_matrix[select_disease_id[row]][select_disease_id[col]]!=0:
        	selected_cos_matrix[row][col]=1
print(selected_cos_matrix.shape)
print(selected_cos_matrix)

np.savetxt('A_disease.txt',selected_cos_matrix,fmt="%.6e")
f_out=open('disease_list_in_A.txt','w')
for i in select_disease:
    f_out.write(i.replace('MESH:','')+'\n')
f_out.close()
