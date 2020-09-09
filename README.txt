INSTALACAO DO StereoisomerIdentifier

O StereoisomerIdentifier é composto por códigos na linguagem C++ e Python.

INSTALACAO C++

Linux
Vá a pasta: 
code/cpp
Digite: 
make

Será gerado um executável de nome: StereoisomerIdentifierRmsd.exe
Mova o executável StereoisomerIdentifierRmsd.exe para a pasta code/python

Em windows a compilação do código C++ pode ser feita com o compilador da sua preferência,
nesse caso, basta adicionar todos os arquivos em sua IDE e compilar.

ATENCAO - o codigo de C++ precisa de uma pasta chamada Stereoisomerlist para rodar, a falta da pasta
pode levar a comportamento imprevisivel.
ESSA PASTA ESTA NO DROPBOX E TEM QUE SER FORNECIDA A PARTE

INSTALACAO PYTHON

Faça do download do anaconda (a internet é rica em tutorias para esse passo).

Crie um ambiente
conda create -n StereoIdent

Ative o ambiente
conda activate StereoIdent

Instale o rdkit
conda install --channel https://conda.anaconda.org/rdkit rdkit

Avisos
- Em uma outra tentativa houve uma mudanca na biblioteca libboost, para corrigir o problema tive que instalar o rdkit com python 3.5.5



EXECUCAO

Ative o ambiente
conda activate StereoIdent

Mova a pasta StereoisomerList para a pasta de trabalho (Essa pasta foi fornecida a parte pelos desenvolvedores)
Na pasta python digite:
python main.py [nome do arquivo]
O resultado sera apresentado no arquivo calculating.csv




ATENCAO - O CODIGO ESTA CONFIGURADO PARA LINUX
- PARA WINDOWS - E NECESSARIO MUDANCAS NO:
               -  StereoiosomerIdentifier.cpp :: filePath (trocar a direcao da barra)
               - CppExecutableHandling.py no windows nao tem esse ponto: ./StereoisomerIdentifierRmsd.exe                           










FORMATO PADRAO DE ARQUIVOS .mol
nome
programa StereoisomerIdentifier
-- blank
natoms nbonds "  0  0  0  0  0  0  0  0999 V2000"
x1 y1 z1 atomlabel1 " 0  0  0  0  0"
.
.
.
xn yn zn atomlabeln " 0  0  0  0  0"
bond1-1 bond1-2 type1 " 0 0 0 " #type always one
.
.
.
bondn-1 bondn-2 typen " 0 0 0 " #type always one
M  END
