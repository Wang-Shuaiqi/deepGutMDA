import re
f=open('../data/MASI_v1.0_download_substanceInfo.txt','r',encoding='utf-8')
f2=open('../data/CTD_chemicals_diseases.csv','r',encoding='utf-8')
f3=open('./disease_list_in_A.txt','r',encoding='utf-8')
f4=open('./drug_disease_interactions.txt','w',encoding='utf-8')
selected_disease=[]
for i in f3:
	i=i.strip('\n')
	selected_disease.append(i)

MASI_substance={}
for i in f:
	i=i.strip('\n').split('\t')
	substance_name=i[1].lower()
	if i[9]!='n.a.':
		MASI_substance[substance_name]=i[9]

useful_substance=set()
MeSH_total={}
for i in f2:
	j= re.split(',(?=([^\"]*\"[^\"]*\")*[^\"]*$)',i)
	if len(j) != 0:
		chem_name=j[0].replace('\'','')
		if chem_name.lower() in  MASI_substance:
			useful_substance.add(chem_name.lower())
			meshid=j[8].replace('MESH:','')
			if meshid in selected_disease:
				if meshid not in MeSH_total:
					MeSH_total[meshid]=set()
					MeSH_total[meshid].add(chem_name.lower())
				else:
					MeSH_total[meshid].add(chem_name.lower())

for disease in MeSH_total:
	for substance in MeSH_total[disease]:
		f4.write(MASI_substance[substance]+'\t'+disease+'\n')
