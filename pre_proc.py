import sys
import re

three2one_dict = {'HSD': 'H', 'HSP': 'H', 'HSE': 'H',  'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}
amino_seq = ''
angles_seq = ''

for l in sys.stdin.readlines()[1:]:
    elems = re.split(r" +", l)
    if elems[0] == '':
        elems.remove('')
    ## print(elems)
    angles_seq += elems[4] + ',' + elems[5]
    amino_seq = amino_seq + three2one_dict[elems[2]]

print(amino_seq)
print()
print(angles_seq)