from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Descriptors import NumRadicalElectrons
from IPython.display import Image, display

def rxn_smarts(reactant_smiles, product_smiles):
    rxn = rdChemReactions.ChemicalReaction()
    for smi in reactant_smiles:
        mol = Chem.MolFromSmiles(smi)
        rxn.AddReactantTemplate(mol)
        
    for smi in product_smiles:
        mol = Chem.MolFromSmiles(smi)
        rxn.AddProductTemplate(mol)

    return rdChemReactions.ReactionToSmarts(rxn)

def run_rxn(sma, reactant_smiles):
    rxn = AllChem.ReactionFromSmarts(sma)
    reactant_mol = tuple( [Chem.AddHs(Chem.MolFromSmiles(smi)) 
                           for smi in reactant_smiles]  )
    products = rxn.RunReactants(reactant_mol)

    possible_products = []
    for prod in products:
        prod_each_mol = []
        for x in prod:
            tmp_mol = Chem.MolFromSmiles(Chem.MolToSmiles(x))
            for atom in tmp_mol.GetAtoms():
                atom.SetNoImplicit(True) # otherwise they just attach H to radical center and remove radical
            smi_tmp = Chem.MolToSmiles(tmp_mol, allHsExplicit = True)
            smi_can = Chem.MolToSmiles(Chem.MolFromSmiles(smi_tmp))
            #print(NumRadicalElectrons(Chem.MolFromSmiles(smi_can)))
            prod_each_mol.append(smi_can) #annoying, but this is the way to preserve a radical center
        possible_products.append(tuple(sorted(prod_each_mol)))   
    possible_products = list(set(possible_products))
    
    return possible_products

def visualize_rxn(rxn_sma, filename):
    rxn = AllChem.ReactionFromSmarts(rxn_sma,useSmiles=True)
    d2d = Draw.MolDraw2DCairo(800,300)
    
    #s = d2d.drawOptions()
    #s.addAtomIndices = True
    
    d2d.DrawReaction(rxn,highlightByReactant=False)
    d2d.FinishDrawing()
    png = d2d.GetDrawingText()
    open(filename,'wb+').write(png)
    return png

def get_rxn_type(rxn_smiles):
    rxn_type = ''
    r_smi, p_smi = rxn_smiles.split('>>')
    
    r_mol = Chem.MolFromSmiles(r_smi)
    r_D = Chem.rdmolops.GetDistanceMatrix(r_mol)
    
    #print(rxn_smiles)
    for atom in r_mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetNumRadicalElectrons() == 1:
            #there is only one carbon whose distance to radical O is two
            R_OO_carbon_index = list(r_D[atom.GetIdx()]).index(2)
            break
    R_OO_carbon_neighbors = r_mol.GetAtomWithIdx(R_OO_carbon_index).GetNeighbors()
    neighbor_symbols = [atom.GetSymbol() for atom in R_OO_carbon_neighbors]
    
    R_OO_carbon_deg = len(R_OO_carbon_neighbors) - 1 # one is OO
    rxn_type += str(R_OO_carbon_deg)
    
    if neighbor_symbols.count('O') >= 2:
        rxn_type += 'a'
    
    rxn_type += '-'
    ######################
    p_mol = Chem.MolFromSmiles(p_smi)
    for atom in p_mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetNumRadicalElectrons() == 1:
            P_Q_carbon_index = atom.GetIdx()
            break
    P_Q_carbon_neighbors = p_mol.GetAtomWithIdx(P_Q_carbon_index).GetNeighbors()
    neighbor_symbols = [atom.GetSymbol() for atom in P_Q_carbon_neighbors]
    
    P_Q_carbon_deg = len(P_Q_carbon_neighbors)
    rxn_type += str(P_Q_carbon_deg)
    
    if neighbor_symbols.count('O') >= 1:
        rxn_type += 'a'

    return rxn_type
