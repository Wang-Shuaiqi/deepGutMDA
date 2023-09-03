f1=open('../data/microbe_in_microphenoDB_id.txt','r')
f2=open('../data/disease_in_microphenoDB_id.txt','r')
f=open('../data/microphenoDB_core_table.txt','r')
f_out=open('microbe_disease_interactions.txt','w')
microbe_name2id={}
disease_name2id={}
for i in f1:
	i=i.strip('\n').split('\t')
	if i[1]!='':
		microbe_name2id[i[0]]=i[1]

for i in f2:
	i=i.strip('\n').split('\t')
	if i[1]!='Error':
		disease_name2id[i[0]]=i[1]

# f_out.write('microbe_id'+'\t'+'disease_id'+'\n')
for i in f:
	i=i.strip('\n').split('\t')
	if i[1] in microbe_name2id and i[2] in disease_name2id:
		f_out.write(microbe_name2id[i[1]]+'\t'+disease_name2id[i[2]]+'\n')

