import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.sparse as sp
import random

import keras
import tensorflow as tf
from keras import Model, optimizers
from keras.callbacks import EarlyStopping
from keras.layers import Input, Dense, Flatten, Dropout
from keras.optimizers import Adam
from keras.regularizers import l2

import spektral
from spektral.layers import GCNConv as GraphConv
from spektral.utils import normalized_laplacian
from keras.utils import np_utils
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score


# gpus = tf.config.experimental.list_physical_devices(device_type='GPU')
# tf.config.experimental.set_visible_devices(devices=gpus[0], device_type='GPU')

drugs = pd.read_table('./drug_list_in_A.txt',header  = None)
microbes = pd.read_table('./microbe_list_in_A.txt',header  = None)
drug_list =[str(i) for i in list(drugs.iloc[:,0])]
microbe_list=[str(i) for i in list(microbes.iloc[:,0])]

A_drug = np.loadtxt('./A_drug.txt')
A_microbe = np.loadtxt('./A_microbe.txt')
laplacian_A_drug = normalized_laplacian(A_drug)
laplacian_A_microbe = normalized_laplacian(A_microbe)
adj=np.vstack((np.hstack((laplacian_A_drug,np.zeros((1041,960)))),np.hstack((np.zeros((960,1041)),laplacian_A_microbe))))
A_drug = np.loadtxt('./drug_features.txt')
A_microbe = np.loadtxt('./microbe_features.txt')

interaction_pair=[]
matrix_drug2microbe=np.zeros((len(drug_list),len(microbe_list)))
matrix_microbe2drug=np.zeros((len(microbe_list),len(drug_list)))
with open('./drug_microbe_interactions.txt','r') as f:
	pair_list=f.read().split('\n')
	pair_list.remove('')
	random.shuffle(pair_list)

for each in pair_list:
    j=each.split('\t')
    microbe=j[1]
    drug=j[0]
    if microbe in microbe_list and drug in drug_list:
        matrix_drug2microbe[drug_list.index(drug)][microbe_list.index(microbe)]=1
        matrix_microbe2drug[microbe_list.index(microbe)][drug_list.index(drug)]=1
        interaction_pair.append(microbe+'_'+drug)

X=[]
Y=[]
Train_data_x=[]
Train_data_y=[]
for md_pair_index in range(len(interaction_pair)):
    microbe, drug  =  interaction_pair[md_pair_index].split('_')
    microbe_microbe=A_microbe[microbe_list.index(microbe)]
    microbe_drug=matrix_microbe2drug[microbe_list.index(microbe)]
    drug_drug=A_drug[drug_list.index(drug)]
    drug_microbe=matrix_drug2microbe[drug_list.index(drug)]
    feature = np.vstack((np.hstack((microbe_drug,microbe_microbe)),np.hstack((drug_drug,drug_microbe)))).T
    X.append(feature)
    Y.append(float(1))
    Train_data_x.append(microbe+'_'+drug)
    Train_data_y.append('1')
    n=1
    while n>0 :
        random_drug = random.sample(drug_list, 1)[0]
        random_microbe = random.sample(microbe_list, 1)[0]
        pair=random_microbe+'_'+random_drug
        if pair not in interaction_pair:
            microbe_microbe=A_microbe[microbe_list.index(random_microbe)]
            microbe_drug=matrix_microbe2drug[microbe_list.index(random_microbe)]
            drug_drug=A_drug[drug_list.index(random_drug)]
            drug_microbe=matrix_drug2microbe[drug_list.index(random_drug)]
            feature = np.vstack((np.hstack((microbe_drug,microbe_microbe)),np.hstack((drug_drug,drug_microbe)))).T
            X.append(feature)
            Y.append(float(0))
            Train_data_x.append(pair)
            Train_data_y.append('0')
            n=n-1
X=np.array(X)
Y=np.array(Y).reshape(len(Y),1)
Train_data_x=np.array(Train_data_x)
Train_data_y=np.array(Train_data_y)



X_score=[]
for i in drug_list:
    for j in microbe_list:
        microbe_microbe=A_microbe[microbe_list.index(j)]
        microbe_drug=matrix_microbe2drug[microbe_list.index(j)]
        drug_drug=A_drug[drug_list.index(i)]
        drug_microbe=matrix_drug2microbe[drug_list.index(i)]
        feature = np.vstack((np.hstack((microbe_drug,microbe_microbe)),np.hstack((drug_drug,drug_microbe)))).T
        X_score.append(feature)
X_score=np.array(X_score)


kf = KFold(n_splits=5,shuffle=False)
n=0
for train_index , test_index in kf.split(X):
    n+=1
    Train_index=train_index[:int(len(train_index)*0.8)]
    Val_index=train_index[int(len(train_index)*0.8):]
    X_train,Y_train=X[Train_index],Y[Train_index]
    X_val,Y_val=X[Val_index],Y[Val_index]
    X_test,Y_test= X[test_index],Y[test_index]
    np.savetxt('./2N_output/x_test%d.txt'%n,Train_data_x[test_index].reshape(test_index.shape[0],1),fmt ='%s')
    np.savetxt('./2N_output/y_test%d.txt'%n,Train_data_y[test_index].reshape(test_index.shape[0],1),fmt ='%s')
    N= X_train.shape[-2]
    F = X_train.shape[-1]
    n_out = Y_train.shape[-1]
    X_in = Input(shape=(N,F))
    graph_conv = GraphConv(64,activation='relu',kernel_regularizer=l2(1e-6),use_bias=True)([X_in, adj])
    graph_conv = GraphConv(64,activation='relu',kernel_regularizer=l2(1e-6),use_bias=True)([graph_conv, adj])
    flatten = Flatten()(graph_conv)
    fc = Dense(128, activation='relu')(flatten)
    fc = Dense(64, activation='relu')(fc)
    fc = Dense(64, activation='relu')(fc)
    dropout=Dropout(0.5)(fc)
    output = Dense(n_out, activation='sigmoid')(dropout)
    model = Model(inputs= X_in, outputs=output)
    model.compile(optimizer=tf.keras.optimizers.Adam(lr=1e-5),loss='binary_crossentropy',metrics=['accuracy'])
    model.summary()
    validation_data = (X_val,Y_val)
    callbacks = [EarlyStopping(monitor='val_accuracy', patience=20)]
    history = model.fit(X_train,Y_train ,batch_size=64,callbacks=callbacks,validation_data=validation_data,epochs=100)
    model.save('./2N_output/model%d.h5'%n)
    Y_predict = model.predict(X_score).reshape(1041,960)
    np.savetxt('./2N_output/scores%d.txt'%n,Y_predict)

