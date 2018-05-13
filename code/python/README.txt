rodando:

Abrir o Anaconda Prompt e deslocar para essa pasta

executar:

conda activate test-rdkit

executar:

python main.py


>>> s = '12abcd405'
>>> result = ''.join([i for i in s if not i.isdigit()])
>>> result
'abcd'



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