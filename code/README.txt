

Instalando python no windows

P�gina de refer�ncia:http://www.rdkit.org/docs/Install.html

Necess�rio instalar o miniconda

Instalacao do RDkit:

1. Abrir o Anaconda
2. conda create --name test-rdkit --channel https://conda.anaconda.org/rdkit rdkit python=3.6.1
3. Ativacao do ambiente:
   - conda activate test-rdkit
4. Instala��o do spyder (editor) dentro do ambiente test-rdkit
   - conda install spyder

Instalacao do matplotlib

conda create --name mpl33 python=3.3 matplotlib ipython-notebook
conda activate mpl33
conda install pandas
conda install scikit-learn

Para rodar o pca.py

conda activate mpl33
ipython pca.py

