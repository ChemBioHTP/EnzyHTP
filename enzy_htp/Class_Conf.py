__doc__ = """
This module make configuration file and deploy input file for different next-stage software. (Amber, Gaussian(QM/ONIOM))
Take all job set up that not relate to the PDB itself
-------------------------------------------------------------------------------------
Config --- Amber
        |- Gaussian
        |- Multiwfn
-------------------------------------------------------------------------------------
"""
from math import exp
import re
from helper import mkdir, line_feed


class Config:
    # >>>>>> Resource <<<<<<
    # -----------------------------
    # Cores available (used in MD(Amber) and QM(Gaussian) calculations)
    #
    n_cores = 24
    # -----------------------------
    # Per core memory in MB
    #
    max_core = 2000
    # -----------------------------
    # debug info level
    # 0: No info. 1: normal debug (warning) 2: verbose debug (running log)
    debug = 1
    # -----------------------------
    # command line name for cpu parallel computing
    #
    PC_cmd = "mpirun"

    @classmethod
    def get_PC_cmd(cls):
        if cls.PC_cmd == "mpirun":
            return "mpirun -np " + str(cls.n_cores)
        return cls.PC_cmd

    # >>>>>> Software <<<<<<


class Layer:
    """
    set / use preset of oniom layer
    ------------
    PDB: related PDB object
    layer: list of layer's list of atom indexes
    """

    def __init__(self, PDB_obj, atom_lists, if_set=0):
        """
        general way to assign layer: list of atom indexes
        """
        self.PDB = PDB_obj
        if not if_set:
            self.layer = []

            if len(atom_lists) not in [2, 3]:
                raise Exception(
                    "Layer.__init__: only support 2 or 3 layers. Input layers: "
                    + str(len(atom_lists))
                )

            for layer in atom_lists:
                atoms = []
                for atom_str in layer.split(","):
                    if "-" in atom_str:
                        if re.match(r"[0-9]+-(?:[0-9]+|L)", atom_str) == None:
                            raise Exception(
                                "Wrong layer syntax of atom indexes. e.g.: 1,2,3-5"
                            )
                        ID = atom_str.split("-")
                        for i, j in enumerate(ID):
                            if j == "L":
                                ID[i] = PDB_obj.get_last_A_id()
                        atoms.extend(range(int(ID[0]), int(ID[1]) + 1))
                    else:
                        atoms.append(int(atom_str))
                self.layer.append(atoms)
        else:
            # for preset
            self.layer = atom_lists

        if Config.debug >= 1:
            print("current layer atoms: ")
            for i, atoms in enumerate(self.layer):
                print(str(i) + ":", len(atoms))
        if Config.debug >= 2:
            print("current layer atoms: " + repr(self.layer))

    @classmethod
    def preset(cls, PDB_obj, set_id, lig_list=[]):
        """
        preset layer settings for ONIOM
        ------------
        set_id =
                (two layers)
                1: Substrate only
                2: Key ligands (specify ligand index or name; all ligand by default)
                3: Substrate and key residues (manually assigned)
                4: Substrate and all residues/ligands within a assigned radius (Need to keep consistant molecular number)
                5: Based on some parameters to select the QM region
        lig_list: (set_id = 2) key ligand index
        """
        layer_atoms = []
        PDB_obj.get_stru()
        stru = PDB_obj.stru

        if set_id not in [1, 2, 3, 4, 5]:
            raise Exception(
                "Only support 1-5 set_id now. You are entering: " + str(set_id)
            )

        if set_id == 1:
            # TODO: connect with the database and recognize the substrate in the future
            pass

        if set_id == 2:
            if lig_list == []:
                h_atoms = []
                lig_names = []
                for lig in stru.ligands:
                    for atom in lig:
                        h_atoms.append(atom.id)
                    lig_names.append(repr((lig.name, lig.id)))
                layer_atoms.append(h_atoms)
                # log
                if Config.debug >= 1:
                    print(
                        "Layer.preset(set_id=2): No id in lig_list, set all ligands to high layer by default"
                    )
                    print(" ".join(lig_names))
            else:
                h_atoms = []
                for lig_inp in lig_list:
                    if type(lig_inp) == int:
                        for lig in stru.ligands:
                            if lig.id == lig_inp:
                                for atom in lig:
                                    h_atoms.append(atom.id)
                    if type(lig_inp) == str:
                        for lig in stru.ligands:
                            if lig.name == lig_inp:
                                for atom in lig:
                                    h_atoms.append(atom.id)

            l_atoms = list(set(stru.get_atom_id()).difference(set(h_atoms)))

            # debug
            # a = stru.get_atom_id()
            # b = set(a)
            # for i in b:
            # 	count = 0
            # 	for j in a:
            # 		if j == i:
            # 			count +=1
            # 	if count > 1:
            # 		print(i, ':', count)

            layer_atoms = [h_atoms, l_atoms]

        if set_id == 3:
            pass

        return cls(PDB_obj, layer_atoms, if_set=1)

    """
	====
	Special Method
	====
	"""

    def __getitem__(self, key: int):
        """
        Chain_obj[int]: Chain_obj.residues[int]
        -----
        use residue index within the chain, start from 0
        """
        return self.layer[key]

    def __len__(self):
        return len(self.layer)
