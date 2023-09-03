f=open('../data/MASI_v1.0_download_substanceInfo.txt','r',encoding='utf-8')
f1=open('../data/MASI_v1.0_download_microbeSubstanceInteractionRecords_ver20200928.txt','r',encoding='utf-8')
f2=open('drug_microbe_interactions.txt','w')
MASIid2pubchemid={}

for i in f:
	i=i.strip('\n').split('\t')
	if i[9]!='n.a.':
		MASIid2pubchemid[i[0]]=i[9]

# f2.write('pubchemid'+'\t'+'ncbi_id'+'\t'+'type'+'\n')
for i in f1:
	i=i.strip('\n').split('\t')
	if i[5] in MASIid2pubchemid and i[3]!='n.a.':
		f2.write(MASIid2pubchemid[i[5]]+'\t'+i[3]+'\t'+i[1]+'\n')

