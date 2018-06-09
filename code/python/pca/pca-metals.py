import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import plotly.plotly as py
import plotly.tools as tls

df = pd.read_csv('metal-statistics.csv',parse_dates=True,sep=r';')
#df.drop(df.index[[46]],inplace=True) #drop if needed

feature = df.columns.values
feature = feature[2:len(feature)]

x = df.loc[:, feature].values
x = StandardScaler().fit_transform(x)
y = df.iloc[:,0].values


pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents,columns = ['principal component 1', 'principal component 2'])
finalDf = pd.concat([principalDf, df.iloc[:,1]], axis = 1)
print(pca.explained_variance_ratio_)

p1 = finalDf.iloc[:,0]
p2 = finalDf.iloc[:,1]

color = df.iloc[:,0]

labels = df.iloc[:,1]

fig = plt.figure(1)
plot = fig.add_subplot(111)
plot.tick_params(axis='both', which='major', labelsize=25)
plot.tick_params(axis='both', which='minor', labelsize=8)
plt.xlabel('PC 1',fontsize=25)
plt.ylabel('PC 2',fontsize=25)

plot.scatter(p1, p2, s=4000, c=color, alpha=0.5)

for i, txt in enumerate(labels):
	if i == 0:#Ag
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.4,p2[i]-0.4),size=40)
	if i == 2:#Au
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i],p2[i]-0.4),size=40)
	if i == 8:#Cu
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.4,p2[i]-0.4),size=40)
	if i == 17:#Hg
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.4,p2[i]-0.4),size=40)
	if i == 21:#La
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.35,p2[i]-0.4),size=40)
	if i == 27:#Ni
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.4,p2[i]-0.4),size=40)
	if i == 29:#Pb
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.4,p2[i]-0.4),size=40)
	if i == 30:#Pd
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-1.0,p2[i]-0.4),size=40)
	if i == 32:#Pt
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.4,p2[i]-0.4),size=40)
	if i == 37:#Sc
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.4,p2[i]-0.4),size=40)
	if i == 45:#U
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.2,p2[i]-0.5),size=40)
	if i == 46:#V
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.2,p2[i]-0.5),size=40)
	if i == 48:#Y
		plot.annotate(txt, xy=(p1[i],p2[i]),xytext=(p1[i]-0.15,p2[i]-0.4),size=40)

plt.show()
