from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit import Chem

'''
#s= 'C1=COC=C1'
s= 'CC[CH2]' #propyl
m = Chem.MolFromSmiles(s)
sma = Chem.MolToSmarts(m, isomericSmiles=True)
sma=sma.replace('#6', 'C').replace('#8', 'O').replace('#15', 'P').replace('#7', 'N')
print(sma)
'''
#rxn = AllChem.ReactionFromSmarts('[C:3][C:2](=[O:1])>>[C:3].[C-:2]#[O+:1]')
#print([Chem.MolToSmiles(x,1) for x in rxn.RunReactants((Chem.MolFromSmiles('CCCC(=O)CC'),))[0]])
#print(Chem.MolToSmiles(Chem.MolFromSmiles('[H]C1CCCCC1[H]')))



rxn = rdChemReactions.ChemicalReaction()

m1 = Chem.MolFromSmiles('COCO[O]')
m1 = Chem.AddHs(m1)
rxn.AddReactantTemplate(m1)
m2= Chem.MolFromSmiles('[C]OCOO')
m2 = Chem.AddHs(m2)
rxn.AddProductTemplate(m2)
print(rdChemReactions.ReactionToSmarts(rxn))




#rxn = AllChem.ReactionFromSmarts('[C:1]-[C:2].[c:3]1:[c:4]:[c:5]:[c:6]:[c:7](-[O:9]):[c:8]:1>>[C:1]=[C:2].[c:3]1:[c:4]:[c:5]:[c:6]:[c:7]:[c:8]:1.[O:9]')

#print([Chem.MolToSmiles(x,1) for x in rxn.RunReactants((Chem.MolFromSmiles('CC'),Chem.MolFromSmiles('c1cccc(O)c1'),))])
#products = rxn.RunReactants((Chem.MolFromSmiles('CCCCCC'),Chem.MolFromSmiles('Cc1cccc(O)c1'),))
#products = rxn.RunReactants((Chem.MolFromSmiles('CC(C)C(C)C'),Chem.MolFromSmiles('Cc1cccc(O)c1')))

#for prod in products:
#    print([Chem.MolToSmiles(x) for x in prod])



'''
# [C:1][O:2][H:3].[H:4][H:5]>>[C:1][H:4].[H:5][O:2][H:3]
#rxn = AllChem.ReactionFromSmarts('[#6:1]~[#8:2].[H:3][H:4]>>[#6:1]([H:3]).[#8:2]([H:4])')
rxn = AllChem.ReactionFromSmarts('[#6:1]1:[#6:2]:[#6:3]:[#6:4]:[#6:5]:[#6:6]:1-[#8:7].[H:8]:[H:9]>>[#6:1]1:[#6:2]:[#6:3]:[#6:4]:[#6:5]:[#6:6]:1.[#8:7]([H:8])[H:9]')
print([Chem.MolToSmiles(x,1) for x in rxn.RunReactants((Chem.MolFromSmiles('c1ccccc1O'),Chem.MolFromSmiles('[HH]'),))])
'''
