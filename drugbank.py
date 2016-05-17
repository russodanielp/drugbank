from rdkit import Chem
import xml.etree.cElementTree as ET

class DrugBankpy:
    """ Initial class to contain all the Drug objects.
    """

    def __init__(self, PATH):
        """ Initialize DrugBank class

        :param PATH: Directory path of XML file
        """
        self.PATH = PATH
        try:
            self.root = ET.parse(self.PATH).getroot()
            self.drugs = list(self.root.findall(".")[0])
        except IOError:
            raise Exception("Incorrect DrugBank XML path specified")

        self.mols = []
        for d in self.drugs:
            tmp = Drug(d)
            if tmp.GetSmiles():
                self.mols.append(Chem.MolFromSmiles(tmp.GetSmiles()))
            elif tmp.GetInChI():
                self.mols.append(Chem.MolFromInchi(tmp.GetInChI()))
        print("Created {0} rdkit mols".format(len(self.mols)))



class Drug:
    """ Individual Drug Object class.  Takes an ElementTree object and initializes and rdkit mol object for each
    drug in DrugBank.
    """
    NAMESPACE = {'drugbank': 'http://www.drugbank.ca'}
    def __init__(self, et_drug):
        """ Initializes a rdkit mol for a DrugBank compound.

        :param et_drug: An ElementTree object of at the parent node of a DrugBank drug.
        """
        self.drug = et_drug

    def GetSmiles(self):
        paths = ["drugbank:experimental-properties/drugbank:property/.[drugbank:kind='SMILES']/drugbank:value",
         "drugbank:calculated-properties/drugbank:property/.[drugbank:kind='SMILES']/drugbank:value"]
        SMILES = None
        for path in paths:
            if self.drug.find(path, Drug.NAMESPACE) is not None:
                SMILES = self.drug.find(path, Drug.NAMESPACE).text
        return SMILES

    def GetInChI(self):
        paths = ["drugbank:experimental-properties/drugbank:property/.[drugbank:kind='InChI']/drugbank:value",
         "drugbank:calculated-properties/drugbank:property/.[drugbank:kind='InChI']/drugbank:value"]
        InChI = None
        for path in paths:
            if self.drug.find(path, Drug.NAMESPACE) is not None:
                InChI = self.drug.find(path, Drug.NAMESPACE).text
        return InChI

