#TODO add documentation
import string
from typing import List, Set
from collections import defaultdict
from biopandas.pdb import PandasPdb

from .atom import Atom
from .chain import Chain
from .residue import Residue
from ..core.logger import _LOGGER
from .ligand import Ligand, residue_to_ligand
from .solvent import Solvent, residue_to_solvent

from ..chemical import one_letters_except, convert_to_one_letter

from .mutate import MutaFlag


class Structure:

    def __init__(self, pandas_pdb: PandasPdb):
        # TODO maybe set the pdb filename as well
        self.chains_ = dict()
        self.residues_ = dict()
        self.pandas_pdb = pandas_pdb
        self.current_pdb = str()
        self.metalatoms = []
        self.ligands = []
        self.solvents = []

    def residue_state(self) -> Set:
        print(self.chains_)
        exit(0)

    def all_possible_mutations(self) -> List[MutaFlag]:
        result = list()
        for cname, chain in self.chains_.items():
            # TODO check that the chain can have this done
            for residue in chain.residues():
                orig_one_letter = convert_to_one_letter(residue.name)
                result.extend(
                    list(
                        map(
                            lambda mut: MutaFlag(orig_residue=orig_one_letter,
                                                 chain_index=cname,
                                                 residue_index=residue.num_,
                                                 target_residue=mut),
                            one_letters_except(orig_one_letter))))
                #print(residue)

        return result

    def get_metalatoms(self):
        """
        get metal from raw chains and clean chains by deleting the metal part
        """
        for cname, chain in self.chains_.items():
            if chain.is_metal():
                self.metalatoms.append(chain)
                del self.chains_[cname]
        if not len(self.metalatoms):
            return
        raise TypeError
        # Break pseudo residues into atoms and convert to Metalatom object
        holders = []
        for pseudo_resi in self.metalatoms:
            for metal in pseudo_resi:
                holders.append(Metalatom.fromAtom(metal))
        metalatoms = holders

        # clean empty chains
        for i in range(len(raw_chains) - 1, -1, -1):
            if len(raw_chains[i]) == 0:
                del raw_chains[i]

        return raw_chains, metalatoms

    def get_ligands(self,
                    ligand_list=None):  # Not sure what type ligand_list is
        """
        get ligand from raw chains and clean chains by deleting the ligand part
        -----
        ligand_list: only record user specified ligand if provided
        Method: Assume metal/ligand/solvent can only be in a seperate chain (or it can not be distinguish from artificial residues.)
                - delete names from rd_non_ligand_list
                + get names from ligand_list only if provided.
        """
        bad_chains = []
        ligands = []
        for cname, chain in self.chains_.items():

            if chain.is_HET():
                continue

            for idx, residue in enumerate(chain.residues()[::-1]):
                if ligand_list is not None and residue.name in ligand_list:
                    _LOGGER.info(
                        f"Structure: Found user assigned ligand in raw {chain.name()} {residue.name} {residue.num_}"
                    )
                    ligands.append(residue)
                    del chain[idx]
                elif not residue.is_rd_solvent(
                ) and not residue.is_rd_non_ligand():
                    _LOGGER.warn(
                        f"Structure: Found ligand in raw {chain.name()} {residue.name} {residue.num_}"
                    )
                    ligands.append(residue)
                    del chain[idx]

            if chain.empty():
                bad_chains.append(cname)

        for bc in bad_chains:
            del self.chains_[bc]

        self.ligands = list(map(residue_to_ligand, ligands))

    def get_solvents(self):
        """
        get solvent from raw chains and clean chains by deleting the solvent part
        -----
        Method: Assume metal/ligand/solvent can anywhere. Base on rd_solvent_list
        """
        solvents = []
        bad_chains = []
        for cname, chain in self.chains_.items():
            for idx, residue in enumerate(chain.residues()[::-1]):
                if residue.is_rd_solvent():
                    _LOGGER.warn(
                        f"Structure: found solvent in raw {residue.name} {residue.id}"
                    )
                    solvents.append(residue)
                    del chain[idx]

            if chain.empty():
                bad_chains.append(cname)
        # Convert pseudo residues to Ligand object
        for bc in bad_chains:
            del self.chains_[bc]

        self.solvents = list(map(residue_to_solvent, solvents))

    def get_metal_center(self):
        """
        Extract metal centers from metalatoms. Judged by the MetalCenter_map
        save to self.metal_centers
        return self.metal_centers
        """
        self.metal_centers = list(
            filter(lambda ma: ma.is_metal_center(), self.metalatoms))
        # TODO maybe this should be copied?
        return self.metal_centers

    def get_art_resi(self):
        """
        find art_resi
        """
        pass

    def add(self, obj, id=None, sort=0):
        _LOGGER.warn(f"Structure.add has not been implemented yet!!!" "")
        return
        """
        1. judge obj type (go into the list)
        2. assign parent
        3. id
        if sort:
            clean original id (use a place holder to represent last)
        if not None:
            assign id
        if sort and not None:
            mark as id+i
        4. add to corresponding list
                                                    sort
                  |     |               0             |                     1               |
         assigned |  0  |   keep (for direct output)  | clean (for sort)                    |
                  |  1  |  assign (for direct output) | mark  (for relative order in sort ) |
        """
        # list
        if type(obj) == list:

            obj_ele = obj[0]

            if (type(obj_ele) != Chain and type(obj_ele) != Metalatom and
                    type(obj_ele) != Ligand and type(obj_ele) != Solvent):
                raise TypeError(
                    "structure.Add() method only take Chain / Metalatom / Ligand / Solvent"
                )

            # add parent and clean id (if sort) assign id (if assigned) leave mark if sort and assigned
            #                         sort
            #          |     |    0     |   1   |
            # assigned |  0  |   keep   | clean |
            #          |  1  |  assign  | mark  |
            for i in obj:
                i.set_parent(self)
                if sort:
                    if id != None:
                        i.id = str(id) + "i"  # str mark
                    else:
                        i.id = id  # None
                else:
                    if id != None:
                        i.id = id

            if type(obj_ele) == Chain:
                self.chains.extend(obj)
            if type(obj_ele) == Metalatom:
                self.metalatoms.extend(obj)
            if type(obj_ele) == Ligand:
                self.ligands.extend(obj)
            if type(obj_ele) == Solvent:
                self.solvents.extend(obj)

        # single building block
        else:
            if (type(obj) != Chain and type(obj) != Metalatom and
                    type(obj) != Ligand and type(obj) != Solvent):
                raise TypeError(
                    "structure.Add() method only take Chain / Metalatom / Ligand / Solvent"
                )

            obj.set_parent(self)
            if sort:
                if id != None:
                    obj.id = str(id) + "i"  # str mark
                else:
                    obj.id = id  # None
            else:
                if id != None:
                    obj.id = id

            if type(obj) == Chain:
                self.chains.append(obj)
            if type(obj) == Metalatom:
                self.metalatoms.append(obj)
            if type(obj) == Ligand:
                self.ligands.append(obj)
            if type(obj) == Solvent:
                self.solvents.append(obj)

        if sort:
            self.sort()

    def sort(self, if_local=0):
        _LOGGER.warn(f"Structure.sort has not been implemented yet!!")
        return
        """
        assign index according to current items
        chain.id
        resi.id
        atom.id
        -----------
        Chain/Residue level: 
            Base on the order of the old obj.id 
            and potential insert mark from add (higher than same number without the mark)
            *if added object has same id and is not assigned with a insert mark -- place after a original one.
        Atom level:
            base on the parent order (parent.id):
            chains -> metalatoms -> ligands
            residue.id within each above.
            list order within each residues.
        """
        if if_local:
            # sort chain order
            self.chains.sort(key=lambda chain: chain.id)
            # rename each chain
            for index, chain in enumerate(self.chains):
                chain.id = chr(65 + index)  # Covert to ABC using ACSII mapping
                # sort each chain
                chain.sort()

            # sort ligand // Do I really need?
            for ligand in self.ligands:
                ligand.sort()  # Do nothing
        else:
            r_id = 0
            a_id = 0
            for chain in self.chains:
                for res in chain:
                    r_id += 1
                    res.id = r_id
                    for atom in res:
                        a_id += 1
                        atom.id = a_id
            for metal in self.metalatoms:
                a_id += 1
                metal.id = a_id
            for lig in self.ligands:
                r_id += 1
                lig.id = r_id
                for atom in lig:
                    a_id += 1
                    atom.id = a_id
            for sol in self.solvents:
                r_id += 1
                sol.id = r_id
                for atom in sol:
                    a_id += 1
                    atom.id = a_id

    def build(self, path, ff="AMBER", forcefield="ff14SB", keep_id=0):
        """
        build PDB after the change based on the chosen format and forcefield
        - line based on atom and contain chain index and residue index
        ----------------------------
        ff = 
        AMBER (standard amber format: from tleap examples)
            - resi and atom indexes start from 1 and DO NOT reset reaching a new chain. 
            - use atom and residue names from amber force field.
            - place metal, ligand, solvent in seperate chains (seperate with TER)
            - ligand -> metal -> solvent order
            * do not sort atomic order in a residue like tleap does.
        """
        with open(path, "w") as of:
            if ff == "AMBER":
                if not keep_id:
                    a_id = 0
                    r_id = 0
                    for chain in self.chains_.values():
                        # write chain
                        for resi in chain.residues():
                            r_id = r_id + 1
                            for atom in resi.atom_list(
                            ):  # TODO fix this horrible name
                                a_id = a_id + 1  # current line index
                                line = atom.build(a_id=a_id,
                                                  r_id=r_id,
                                                  ff=ff,
                                                  forcefield=forcefield)
                                of.write(line + '\n')
                        # write TER after each chain
                        of.write("TER\n")

                    c_id = chr(len(self.chains_) + 64)

                    for ligand in self.ligands:
                        r_id = r_id + 1
                        c_id = chr(ord(c_id) + 1)

                        for atom in ligand:
                            a_id = a_id + 1
                            line = atom.build(
                                a_id=a_id,
                                r_id=r_id,
                                c_id=c_id,
                                ff=ff,
                                forcefield=forcefield,
                            )
                            of.write(line + '\n')
                        of.write("TER" + line_feed)

                    for metal in self.metalatoms:
                        a_id = a_id + 1
                        r_id = r_id + 1
                        c_id = chr(ord(c_id) + 1)

                        line = metal.build(
                            a_id=a_id,
                            r_id=r_id,
                            c_id=c_id,
                            ff=ff,
                            forcefield=forcefield,
                        )
                        of.write(line + '\n')
                        of.write("TER" + line_feed)

                    if len(self.solvents) != 0:
                        c_id = chr(ord(c_id) +
                                   1)  # same chain_id for all solvent
                        for solvent in self.solvents:
                            r_id = r_id + 1
                            for atom in solvent:
                                a_id = a_id + 1
                                line = atom.build(
                                    a_id=a_id,
                                    r_id=r_id,
                                    c_id=c_id,
                                    ff=ff,
                                    forcefield=forcefield,
                                )
                                of.write(line + '\n')
                        of.write("TER" + line_feed)

                else:
                    for chain in self.chains:
                        # write chain
                        for resi in chain:
                            for atom in resi:
                                line = atom.build(ff=ff, forcefield=forcefield)
                                of.write(line + '\n')
                        # write TER after each chain
                        of.write("TER" + line_feed)

                    c_id = chr(len(self.chains) + 64)

                    for ligand in self.ligands:
                        c_id = chr(ord(c_id) + 1)
                        for atom in ligand:
                            line = atom.build(c_id=c_id,
                                              ff=ff,
                                              forcefield=forcefield)
                            of.write(line + '\n')
                        of.write("TER" + line_feed)

                    for metal in self.metalatoms:
                        c_id = chr(ord(c_id) + 1)
                        line = metal.build(c_id=c_id,
                                           ff=ff,
                                           forcefield=forcefield)
                        of.write(line + '\n')
                        of.write("TER" + line_feed)

                    if len(self.solvents) != 0:
                        c_id = chr(ord(c_id) + 1)  # chain_id for all solvent
                        for solvent in self.solvents:
                            for atom in solvent:
                                line = atom.build(c_id=c_id,
                                                  ff=ff,
                                                  forcefield=forcefield)
                                of.write(line + '\n')
                        of.write("TER" + line_feed)

            if ff == "XXX":
                # place holder
                pass

            of.write("END\n")

    def build_ligands(self,
                      dir,
                      ft="PDB",
                      ifcharge=0,
                      c_method="PYBEL",
                      ph=7.0,
                      ifname=0,
                      ifunique=0):
        """
        build files for every ligand in self.ligands
        -------
        dir      : output dir. (e.g. File path for ligand i is $dir/ligand_i.pdb)
        ft       : file type / now support: PDB(default)
        ifcharge : if calculate net charge info. (do not before add H)
        c_method : method determining the net charge (default: PYBEL)
        ph       : pH value used for determine the net charge
        ifname   : export residue name if 1 (default: 0)
        ifunique : 1: only build one ligand if there's multiple same ones. 0: build every ligand
        """
        out_ligs = []
        # TODO ph check
        l_id = 0
        lig_list = self.ligands
        # Only build ligand with same name once
        if ifunique:
            lig_list_uni = {}
            for lig in lig_list:
                lig_list_uni[lig.name] = lig
            lig_list = lig_list_uni.values()

        for lig in lig_list:
            l_id = l_id + 1  # current ligand id
            # make output path
            if dir[-1] == "/":
                dir = dir[:-1]
            if ifunique:
                out_path = dir + "/ligand_" + lig.name + ".pdb"
            else:
                out_path = dir + "/ligand_" + str(
                    l_id) + "_" + lig.name + ".pdb"
            # write
            lig.build(out_path)
            # net charge
            net_charge = None
            if ifcharge:
                if lig.net_charge != None:
                    net_charge = lig.net_charge
                else:
                    net_charge = lig.get_net_charge(method=c_method,
                                                    ph=ph,
                                                    o_dir=dir)

            # record
            if ifname:
                out_ligs.append((out_path, net_charge, lig.name))
            else:
                out_ligs.append((out_path, net_charge))

        return out_ligs

    def build_protein(self, dir, ft="PDB"):
        """
        build only protein and output under the dir
        -------
        dir: out put dir ($dir/protein.pdb)
        ft : file type / now support: PDB(default) 
        """
        # make path
        if dir[-1] == "/":
            dir = dir[:-1]
        out_path = dir + "/protein.pdb"

        # write
        if ft == "PDB":
            with open(out_path, "w") as of:
                a_id = 0
                r_id = 0
                for chain in self.chains:
                    # write chain
                    for resi in chain:
                        r_id = r_id + 1
                        for atom in resi:
                            a_id = a_id + 1  # current line index
                            line = atom.build(a_id=a_id, r_id=r_id)
                            of.write(line)
                    # write TER after each chain
                    of.write("TER" + line_feed)
                of.write("END" + line_feed)
        else:
            raise Exception("Support only PDB output now.")

        return out_path

    def build_metalcenters(self, dir, ft="PDB"):
        """
        build metalcenters only. Use for MCPB parameterization. Deal with donor residue with different protonation states.        ----------
        TODO
        """
        out_paths = []
        return out_paths

    def get_connect(self, metal_fix=1, ligand_fix=1, prepi_path=None):
        """
        get connectivity
        -----------------
        TREATMENT
        chain: based on connectivity map of each atom in each residue
        metalatom:  fix1: treat as isolated atom
                    fix2: connect to donor atom (MCPB?)
        ligand: fix1: use antechamber generated prepin file to get connectivity.
                      according to https://ambermd.org/doc/prep.html the coordniate line will always start at the 11th line after 3 DUMM.
        """
        # san check
        if ligand_fix == 1 and prepi_path == None:
            raise Exception("Ligand fix 1 requires prepin_path.")
        # chain part
        for chain in self.chains:
            for res in chain:
                for atom in res:
                    atom.get_connect()
        for sol in self.solvents:
            for atom in sol:
                atom.get_connect()
        # metal
        for metal in self.metalatoms:
            metal.connect = []
        if metal_fix == 1:
            pass
        if metal_fix == 2:
            raise Exception("TODO: Still working on 2 right now")

        # ligand
        # init
        for lig in self.ligands:
            for atom in lig:
                atom.connect = []
        # fix 1
        if ligand_fix == 1:
            for lig in self.ligands:
                # read prepin for each ligand
                with open(prepi_path[lig.name]) as f:
                    line_id = 0
                    if_loop = 0
                    for line in f:
                        line_id += 1
                        if line.strip() == "":
                            if if_loop == 1:
                                # switch off loop and break if first blank after LOOP encountered
                                if_loop = 0
                                break
                            continue
                        if if_loop:
                            lp = line.strip().split()
                            lig._find_atom_name(lp[0]).connect.append(
                                lig._find_atom_name(lp[1]))
                            continue
                        # loop connect starts at LOOP
                        if line.strip() == "LOOP":
                            if_loop = 1
                            continue
                        # coord starts at 11th
                        if line_id >= 11:
                            lp = line.strip().split()
                            atom_id = int(lp[0]) - 3
                            atom_cnt = int(lp[4]) - 3
                            if atom_cnt != 0:
                                lig[atom_id - 1].connect.append(lig[atom_cnt -
                                                                    1])
                                lig[atom_cnt - 1].connect.append(lig[atom_id -
                                                                     1])

    def get_connectivty_table(self,
                              ff="GAUSSIAN",
                              metal_fix=1,
                              ligand_fix=1,
                              prepi_path=None):
        """
        get connectivity table with atom index based on 'ff' settings:
        ff = GAUSSIAN  -- continuous atom index start from 1, do not seperate by chain
        -------------------
        TREATMENT
        Use original atom.id.
            chain: based on connectivity map of each atom in each residue
            metalatom:  fix1: treat as isolated atom
                        fix2: connect to donor atom (MCPB?)
            ligand: fix1: use antechamber generated prepin file to get connectivity.
        Use 1.0 for all connection.
            Amber force field in gaussian do not account in bond order. (Only UFF does.)
            Note that bond order less than 0.1 do not count in MM but only in opt redundant coordinate.
        """
        connectivty_table = ""
        # get connect for every atom in stru
        self.get_connect(metal_fix, ligand_fix, prepi_path)

        # write str in order
        # Note: Only write the connected atom with larger id
        a_id = 0
        for chain in self.chains:
            for res in chain:
                for atom in res:
                    a_id += 1
                    cnt_line = " " + str(atom.id)
                    # san check
                    if atom.id != a_id:
                        raise Exception("atom id error.")
                    for cnt_atom in atom.connect:
                        if cnt_atom.id > atom.id:
                            cnt_line += " " + str(cnt_atom.id) + " " + "1.0"
                    connectivty_table += cnt_line + line_feed

        for lig in self.ligands:
            for atom in lig:
                a_id += 1
                cnt_line = " " + str(atom.id)
                # san check
                if atom.id != a_id:
                    raise Exception("atom id error.")
                for cnt_atom in atom.connect:
                    if cnt_atom.id > atom.id:
                        cnt_line += " " + str(cnt_atom.id) + " " + "1.0"
                connectivty_table += cnt_line + line_feed

        for atom in self.metalatoms:
            a_id += 1
            cnt_line = " " + str(atom.id)
            # san check
            if atom.id != a_id:
                raise Exception("atom id error.")
            for cnt_atom in atom.connect:
                if cnt_atom.id > atom.id:
                    cnt_line += " " + str(cnt_atom.id) + " " + "1.0"
            connectivty_table += cnt_line + line_feed

        for sol in self.solvents:
            for atom in sol:
                a_id += 1
                cnt_line = " " + str(atom.id)
                # san check
                if atom.id != a_id:
                    raise Exception("atom id error.")
                for cnt_atom in atom.connect:
                    if cnt_atom.id > atom.id:
                        cnt_line += " " + str(cnt_atom.id) + " " + "1.0"
                connectivty_table += cnt_line + line_feed

        return connectivty_table

    def protonation_metal_fix(self, Fix):
        """
        return a bool: if there's any metal center
        """
        # try once if not exist
        if self.metal_centers == []:
            self.get_metal_center()
        if self.metal_centers == []:
            print("No metal center is found. Exit Fix.")
            return False

        # start fix
        # get donor atoms and residues
        for metal in self.metal_centers:
            metal.get_donor_residue(method="INC")

            if Fix == 1:
                metal._metal_fix_1()

            if Fix == 2:
                metal._metal_fix_2()

            if Fix == 3:
                metal._metal_fix_3()
        return True

    def get_all_protein_atom(self):
        """
        get a list of all protein atoms
        return all_P_atoms 
        """
        all_P_atoms = []
        for chain in self.chains:
            for residue in chain:
                all_P_atoms.extend(residue.atoms)
        return all_P_atoms

    def get_all_residue_unit(self, ifsolvent=0):
        all_r_list = []
        for chain in self.chains:
            for resi in chain:
                all_r_list.append(resi)

        for resi in self.ligands:
            all_r_list.append(resi)

        for resi in self.metalatoms:
            all_r_list.append(resi)

        if ifsolvent:
            for resi in self.solvents:
                all_r_list.append(resi)

        return all_r_list

    def get_residue(self, id):
        """
        re-search residue id with all residues count togethor from 1.
        ----------
        return a residue object
        """
        all_resi = self.get_all_residue_unit()
        for resi in all_resi:
            if resi.id == int(id):
                return resi

    def get_atom_id(self):
        """
        return a list of id of all atoms in the structure
        """
        atom_id_list = []
        for chain in self.chains:
            for res in chain:
                for atom in res:
                    if atom.id == None:
                        raise Exception("Detected None in chain " +
                                        str(chain.id) + str(res.id) + " " +
                                        atom.name)
                    atom_id_list.append(atom.id)

        for metal in self.metalatoms:
            if metal.id == None:
                raise Exception("Detected None in metal " + metal.name)
            atom_id_list.append(metal.id)

        for lig in self.ligands:
            for atom in lig:
                if atom.id == None:
                    raise Exception("Detected None in ligands", res.id,
                                    atom.name)
                atom_id_list.append(atom.id)

        for sol in self.solvents:
            for atom in sol:
                if atom.id == None:
                    raise Exception("Detected None in solvent", res.id,
                                    atom.name)
                atom_id_list.append(atom.id)

        return atom_id_list

    def get_atom_charge(self, prmtop_path):
        """
        requires generate the stru using !SAME! PDB as one that generate the prmtop. 
        """
        pass

    def get_atom_type(self, prmtop_path):
        """
        requires generate the stru using !SAME! PDB as one that generate the prmtop. 
        """
        # get type list
        with open(prmtop_path) as f:
            type_list = []
            line_index = 0

            for line in f:

                line_index = line_index + 1  # current line

                if line.strip() == r"%FLAG POINTERS":
                    format_flag = line_index
                if line.strip() == r"%FLAG AMBER_ATOM_TYPE":
                    type_flag = line_index

                if "format_flag" in dir():
                    if line_index == format_flag + 2:
                        N_atom = int(line.split()[0])
                        del format_flag

                if "type_flag" in dir():
                    if (line_index >= type_flag + 2 and
                            line_index <= type_flag + 1 + ceil(N_atom / 5)):
                        for i in line.strip().split():
                            type_list.append(i)
        # assign type to atom
        for chain in self.chains:
            for res in chain:
                for atom in res:
                    atom.type = type_list[atom.id - 1]
        for atom in self.metalatoms:
            if atom.id == None:
                raise Exception("Detected None in metal " + atom.name)
            atom.type = type_list[atom.id - 1]
        for lig in self.ligands:
            for atom in lig:
                if atom.id == None:
                    raise Exception("Detected None in ligands", res.id,
                                    atom.name)
                atom.type = type_list[atom.id - 1]
        for sol in self.solvents:
            for atom in sol:
                if atom.id == None:
                    raise Exception("Detected None in solvent", res.id,
                                    atom.name)
                atom.type = type_list[atom.id - 1]

    def get_sele_list(self, atom_mask, fix_end="H", prepi_path=None):
        """
        interface with class ONIOM_Frame. Generate a list for sele build. Make sure use same pdb as the one generate the frame.
        ------------
        resi_list: selected residue list
        atom_mask: atom selection with the standard grammer of Amber (incomplete)
        fix_end: fix valence of the cut bond. (default: H)
                - H: add H to where the original connecting atom is.
                    special fix for classical case:
                    "sele by residue" (cut N-C) -- adjust dihedral for added H on N.
                = Interface with write_sele_lines:
                    add {fix_flag+element_mark: coord} in sele_lines 
        ------------
        return a sele list:
        - Fixing atoms are labeled as qm_atom_id-qm_atom_cnt_id-distance
        - backbone atoms are marked as b at the end (for qmcluster charge calculation)
        - other atoms use _ as place holder
        return a sele map:
        (PDB atom id -> QM atom id)
        """
        sele_lines = {}
        # decode atom_mask (maybe in helper later) TODO
        resi_list = atom_mask[1:].strip().split(",")
        all_resi_list = self.get_all_residue_unit()

        # decode and get obj
        sele_stru_objs = []
        for resi in resi_list:
            chain_id = re.match("[A-Z]", resi)
            resi_id = int(re.match("[0-9]+", resi).group(0))
            if chain_id == None:
                for resi in all_resi_list:
                    if resi_id == resi.id:
                        resi_obj = resi
            else:
                chain_id = chain_id.group(0)
                resi_obj = self.chains[int(chain_id) -
                                       65]._find_resi_id(resi_id)

            sele_stru_objs.append(resi_obj)

        # combine the sele
        sele_atoms = []
        for obj in sele_stru_objs:
            for atom in obj:
                sele_atoms.append(atom)

        if fix_end != None:
            self.get_connect(prepi_path=prepi_path)

        # operate on the sele objs
        for atom in sele_atoms:
            # add current atom
            atom.get_ele()
            if type(atom.parent) != Ligand:
                if atom.name in ["C", "CA", "O", "N", "H", "HA"]:
                    sele_lines[str(atom.id) + "b"] = atom.ele
                else:
                    sele_lines[str(atom.id) + "_"] = atom.ele
            else:
                sele_lines[str(atom.id) + "_"] = atom.ele

            if fix_end != None:
                # search for cut bond
                for cnt_atom in atom.connect:
                    if not cnt_atom in sele_atoms:
                        if fix_end == "H":
                            d_XH = X_H_bond_length[atom.name]
                            label = "-".join(
                                (str(atom.id), str(cnt_atom.id), str(d_XH)))
                            fix_atom = "H"
                        if fix_end == "Me":
                            # TODO
                            pass
                        # write to sele_lines
                        sele_lines[label] = fix_atom
        # make sele_map (PDB atom id -> QM atom id)
        sele_map = {}
        for i, key in enumerate(sele_lines.keys()):
            if key[-1] not in "1234567890":
                key = key[:-1]
            sele_map[key] = i + 1

        if Config.debug >= 1:
            print("Selected QM cluster atoms: ")
            print(sele_lines)

        return sele_lines, sele_map

    def get_resi_dist(self, r1, r2, method="mass_center"):
        """
        r1: residue 1. (residue_obj)
        r2: residue 2. (residue_obj)
        method: The method to narrow the residue down to a point.
                - mass_center
                - ...
        """
        if method == "mass_center":
            p1 = r1.get_mass_center()
            p2 = r2.get_mass_center()

        D = get_distance(p1, p2)

        return D

    """
    ====
    Special Method
    ====
    """

    def __len__(self):
        """
        len(obj) = len(obj.child_list)
        """
        return (len(self.chains) + len(self.metalatoms) + len(self.ligands) +
                len(self.solvents))

    def __getitem__(self, i):
        """
        pop four elements with fixed order
        -----------------
        capable with future update
        """

        if i == 0:
            return self.chains
        if i == 1:
            return self.ligands
        if i == 2:
            return self.metalatoms
        if i == 3:
            return self.solvents
        if i > 3:
            raise StopIteration

    def name_chains(self, mapper: defaultdict) -> None:

        def legal_chain_names(mapper) -> List[str]:
            result = list(string.ascii_uppercase)
            taken = set(list(mapper.keys()))
            result = list(filter(lambda s: s not in taken, result))
            return list(reversed(result))

        key_names = set(list(map(lambda kk: kk.strip(), mapper.keys())))
        if '' not in key_names:
            return mapper
        unnamed = list(mapper[''])
        del mapper['']

        names = legal_chain_names(mapper)
        unnamed = sorted(unnamed, key=lambda r: -r.min_line())
        new_chain: List[Residue] = [unnamed.pop()]

        while len(unnamed):
            new_res = unnamed.pop()
            if new_chain[-1].neighbors(new_res):
                new_chain.append(new_res)
                continue

            new_name = names.pop()
            mapper[new_name] = new_chain
            new_chain = [new_res]

        if len(new_chain):
            new_name = names.pop()
            mapper[new_name] = new_chain

        return mapper

    def build_chains(self) -> None:
        """"""
        # TODO check if there are chain_ids
        self.build_residues()

        mapper = defaultdict(list)
        for res in self.residues_.values():
            mapper[res.chain()].append(res)

        mapper = self.name_chains(mapper)
        # ok this is where we handle missing chain ids
        for chain_name, residues in mapper.items():
            self.chains_[chain_name] = Chain(
                chain_name, sorted(residues, key=lambda r: r.num()))

    def build_residues(self) -> None:
        mapper = defaultdict(list)
        for i, row in self.pandas_pdb.df["ATOM"].iterrows():
            aa = Atom(**row)
            mapper[aa.residue_key()].append(aa)

        for res_key, atoms in mapper.items():
            self.residues_[res_key] = Residue(residue_key=res_key,
                                              atoms=sorted(
                                                  atoms,
                                                  key=lambda a: a.atom_number))


def structure_from_pdb(fname: str) -> Structure:
    """
    extract the structure from PDB path. Capable with raw experimental and Amber format
    ---------
    input = path (or file or file_str)
        split the file_str to chain and init each chain
    ligand_list: ['NAME',...]
        User specific ligand names. Only extract these if provided. 
    ---------
    Target:
    - structure(w/name)  - chain - residue - atom
                        |- metalatom(atom)
                        |- ligand(residue)
                        |- solvent(residue)
    - ... (add upon usage)
    """
    parser = PandasPdb()
    parser.read_pdb(fname)
    result = Structure(parser)
    result.build_chains()
    result.get_metalatoms()
    result.get_ligands()
    result.get_solvents()
    return result
