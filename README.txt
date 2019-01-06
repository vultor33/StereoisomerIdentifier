
INSTALACAO DO StereoisomerIdentifier

O StereoisomerIdentifier � composto por c�digos na linguagem C++ e Python.

INSTALACAO C++

Linux
V� a pasta: 
code/cpp
Digite: 
make

Ser� gerado um execut�vel de nome: StereoisomerIdentifierRmsd.exe
Mova o execut�vel StereoisomerIdentifierRmsd.exe para a pasta code/python

Em windows a compila��o do c�digo C++ pode ser feita com o compilador da sua prefer�ncia,
nesse caso, basta adicionar todos os arquivos em sua IDE e compilar.

INSTALACAO PYTHON

Fa�a do download do anaconda (a internet � rica em tutorias para esse passo).

Crie um ambiente
conda create -n StereoIdent

Ative o ambiente
conda activate StereoIdent

Instale o rdkit
conda install --channel https://conda.anaconda.org/rdkit rdkit


EXECUCAO

Mova a pasta StereoisomerList para a pasta de trabalho

(Essa pasta foi fornecida a parte pelos desenvolvedores)



INSTALACAO DO PYTHON (Windows)

Resumo: Instalar o RDKIT no windows.

P�gina de refer�ncia: http://www.rdkit.org/docs/Install.html
Primeiro eu instalei miniconda.
Abri o conda e executei o seguinte comando:
conda create --name test-rdkit --channel https://conda.anaconda.org/rdkit rdkit python=3.6.1


EXECUTANDO O IDENTIFICADOR

1. Compilar o c�digo na pasta cpp
2. Mover o "StereoisomerIdentifierRmsd.exe" para a pasta do python
3. Abrir o Anaconda Prompt e deslocar para essa pasta
4. Executar:
	conda activate test-rdkit
5. Em seguida, executar:
	python main.py

Ele est� configurado para identificar todos os arquivos .mol2 que estiverem
dentro da pasta.










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
