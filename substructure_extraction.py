import glob
import sys
from rdkit import Chem
from rdkit.Chem import rdFMCS



# Find all .mol files in a current directory

mol_files = glob.glob('*.mol')
mol_list= []

# create a list of input molecules

for mol in mol_files:
    single = Chem.MolFromMolFile(mol)
    mol_list.append(single)

# find maximum common substructure (bond order is set to be flexible)

res = rdFMCS.FindMCS(mol_list, bondCompare=rdFMCS.BondCompare.CompareAny).smartsString
pattern = Chem.MolFromSmarts(res)


# read a molecule in mol format

m = Chem.MolFromMolFile(sys.argv[1])
conf = m.GetConformer()
sub = m.GetSubstructMatch(pattern)

# get coordinates of a maximum common substructure of a molecule

print(len(sub))
print('')
for s in sub:
    coordinates=conf.GetAtomPosition(s)
    print(str(m.GetAtoms()[s].GetSymbol()) + " " + str(coordinates.x) + " " + str(coordinates.y) + " " + str(coordinates.z)
