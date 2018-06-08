import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

df = pd.read_csv('metal-statistics.csv',parse_dates=True,sep=r';')

feature = df.columns.values
feature = feature[2:len(feature)]

x = df.loc[:, feature].values
x = StandardScaler().fit_transform(x)
y = df.iloc[:,0].values

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents,columns = ['principal component 1', 'principal component 2'])
finalDf = pd.concat([principalDf, df.iloc[:,1]], axis = 1)

#print(pca.explained_variance_ratio_)

p1 = finalDf.iloc[:,0]
p2 = finalDf.iloc[:,1]

color = df.iloc[:,0]

labels = df.iloc[:,1]

plt.scatter(p1, p2, s=2000, c=color)

for i, txt in enumerate(labels):
    plt.annotate(txt, (p1[i],p2[i]))

plt.show()
