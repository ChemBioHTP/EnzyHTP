class Residue:

    def __init__(self, residue_key, atoms):
        self.atoms = atoms
        self.residue_key = residue_key
        (chain, name, num) = self.residue_key.split('.')
        self.chain_ = chain
        self.name = name
        self.num_ = int(num)
        # TODO what are the checks that we should be doing here?

    def chain(self):
        return self.chain_

    def num(self):
        return self.num_

    #TODO add operator overloading to move the chain?
