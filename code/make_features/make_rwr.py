from pyrwr.rwr import RWR
import numpy as np
X = np.loadtxt('../drug_similarity.txt', dtype=float, comments='#')
n=X.shape[0]
X_out=np.zeros((n*n,3))
line=0
for i in range(n):
	for j in range(n):
		X_out[line][0]=i
		X_out[line][1]=j
		if X[i][j]>=0:
			X_out[line][2]=X[i][j]
		if X[i][j]<0:
			X_out[line][2]=X[i][j]*-1
		line+=1
np.savetxt('X_drug.txt', X_out)

rwr = RWR()
for each in range(n):
	rwr.read_graph('X_drug.txt', 'directed')
	r = rwr.compute(each)
	if each == 0:
		X_rwr=r
	else:
		X_rwr=np.vstack((X_rwr,r))
np.savetxt('drug_features.txt', X_rwr)
#####
X = np.loadtxt('../microbe_similarity.txt', dtype=float, comments='#')
n=X.shape[0]
X_out=np.zeros((n*n,3))
line=0
for i in range(n):
	for j in range(n):
		X_out[line][0]=i
		X_out[line][1]=j
		if X[i][j]>=0:
			X_out[line][2]=X[i][j]
		if X[i][j]<0:
			X_out[line][2]=X[i][j]*-1
		line+=1
np.savetxt('X_microbe.txt', X_out)

rwr = RWR()
for each in range(n):
	rwr.read_graph('X_microbe.txt', 'directed')
	r = rwr.compute(each)
	if each == 0:
		X_rwr=r
	else:
		X_rwr=np.vstack((X_rwr,r))
np.savetxt('microbe_features.txt', X_rwr)

