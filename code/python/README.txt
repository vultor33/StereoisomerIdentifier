rodando:

Abrir o Anaconda Prompt e deslocar para essa pasta

executar:

conda activate test-rdkit

executar:

python main.py


para rodar precisa do executavel do cpp

e da pasta StereoisomerList - que esta no dropbox : 
!PROJETOS/StereoisomerIdentifier




Formato .mol

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
