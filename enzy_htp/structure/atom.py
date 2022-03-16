import numpy as np


class Atom:

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)
        # self.name = kwargs.get('atom_name',str())
        ##self.x = kwargs.get('x_coord',float())
        ##self.y = kwargs.get('y_coord',float())
        ##self.z = kwargs.get('z_coord',float())
        # self.atom_number = kwargs.get('atom_number', int())
        # self.residue_number = kwargs.get('residue_number', int())
        # self.residue_name = kwargs.get('residue_name', int())
        # self.chain_id = kwargs.get('chain_id', str())

    def __add__(self, offset):
        return Atom(
            name=self.name,
            x=self.x_coord + other,
            y=self.y_coord + other,
            z=self.z_coord + other,
        )

    def __sub__(self, offset):
        return Atom(
            x_coord=self.x - other,
            y_coord=self.y - other,
            z_coord=self.z - other,
            **self.__dict__,
        )

    def residue_key(self):
        return f"{self.chain_id}.{self.residue_name}.{self.residue_number}"

    def distance(self, other):
        return np.sqrt((self.x_coord - other.x_coord)**2 +
                       (self.y_coord - other.y_coord)**2 +
                       (self.z_coord - other.z_coord)**2)

    def __getitem__(self, key):
        return self.__dict__[key]

#class Atom(Child):
#    """
#    -------------
#    initilize from
#    PDB:        Atom.fromPDB(atom_input, input_type='PDB_line' or 'line_str' or 'file' or 'path')
#    raw data:   Atom(atom_name, coord)
#    -------------
#    id
#    name
#    coord = [x, y, z]
#    ele (obtain by method)
#
#    parent # resi
#    -------------
#    Method
#    -------------
#    set_parent
#    Add_id # use only after the whole structure is constructed
#
#    get_connect
#    get_protons
#    get_lp_infos
#    get_bond_end_atom
#    set_byDihedral
#    set_byAngle
#    set_byBond
#
#    gen_line
#
#    get_around
#    get_ele
#    -------------
#    """
#
#    def __init__(self, atom_name: str, coord: list, ff: str, atom_id=None, parent=None):
#        """
#        Common part of init methods: direct from data objects
#        No parent by default. Add parent by action
#        """
#        # set parent to None
#        Child.__init__(self)
#        # add parent if provided
#        if parent != None:
#            self.set_parent(parent)
#        # get data
#        self.name = atom_name
#        self.coord = coord
#        self.id = atom_id
#        self.ff = ff
#
#    @classmethod
#    def fromPDB(cls, atom_input, input_type="PDB_line", ff="Amber"):
#        """
#        generate atom from PDB. Require 'ATOM' and 'HETATM' lines.
#        ---------
#        atom_input = PDB_line (or line_str or file or path)
#        """
#        # adapt general input // converge to a PDB_line (atom_line)
#        if input_type == "path":
#            f = open(atom_input)
#            atom_line = PDB_line(f.read())
#            f.close()
#        if input_type == "file":
#            atom_line = PDB_line(atom_input.read())
#        if input_type == "line_str":
#            atom_line = PDB_line(atom_input)
#        if input_type == "PDB_line":
#            atom_line = atom_input
#
#        # get data
#        atom_name = atom_line.atom_name
#        coord = [atom_line.atom_x, atom_line.atom_y, atom_line.atom_z]
#        atom_id = atom_line.atom_id
#
#        return cls(atom_name, coord, ff, atom_id=atom_id)
#
#    """
#    ====
#    Method
#    ====
#    """
#
#    def get_connect(self):
#        """
#        find connect atom base on:
#        1. AmberMaps.resi_cnt_map
#        2. parent residue name
#        ------------
#        * require standard Amber format (atom name and C/N terminal name)
#        save found list of Atom object to self.connect
#        """
#        self.connect = []
#
#        if self.resi.name in rd_solvent_list:
#            name_list = resi_cnt_map[self.resi.name][self.name]
#        else:
#            r = self.resi
#            r1 = self.resi.chain.residues[0]
#            rm1 = self.resi.chain.residues[-1]
#
#            if r == r1:
#                # N terminal
#                name_list = resi_nt_cnt_map[self.resi.name][self.name]
#            else:
#                if r == rm1:
#                    # C terminal
#                    name_list = resi_ct_cnt_map[self.resi.name][self.name]
#                else:
#                    name_list = resi_cnt_map[self.resi.name][self.name]
#
#        for name in name_list:
#            try:
#                if name not in ["-1C", "+1N"]:
#                    cnt_atom = self.resi._find_atom_name(name)
#                if name == "-1C":
#                    cnt_resi = self.resi.chain._find_resi_id(self.resi.id - 1)
#                    cnt_atom = cnt_resi._find_atom_name("C")
#                if name == "+1N":
#                    cnt_resi = self.resi.chain._find_resi_id(self.resi.id + 1)
#                    cnt_atom = cnt_resi._find_atom_name("N")
#                self.connect.append(cnt_atom)
#            except IndexError:
#                if Config.debug >= 1:
#                    print(
#                        "WARNING: "
#                        + self.resi.name
#                        + str(self.resi.id)
#                        + " should have atom: "
#                        + name
#                    )
#                else:
#                    pass
#
#    def get_type(self):
#        if self.resi.name in rd_solvent_list:
#            self.type = G16_label_map[self.parent.name][self.name]
#        else:
#            r = self.resi
#            r1 = self.resi.chain.residues[0]
#            rm1 = self.resi.chain.residues[-1]
#            if r == r1 and self.name in ["H1", "H2", "H3"]:
#                self.type = "H"
#            else:
#                if r == rm1 and self.name == "OXT":
#                    self.type = "O2"
#                else:
#                    self.type = G16_label_map[self.parent.name][self.name]
#        return self.type
#
#    def get_pseudo_H_type(self, old_atom, method=1):
#        """
#        determine the connecting pseudo H atom type base on the environment for ONIOM input
#        A low layer atom will be replaced to a pesudo H for "model-low" layer modeling
#        self    : the atom in the "model" layer
#        old_atom: the original connecting atom replacing by H.
#        method  :
#        1 ----- : Only consider some main cases:
#                    - truncating CA-CB
#                    - truncating C-N
#                    (abort for other cases)
#        -----------
#        TARGET
#        -----------
#        1. provide reasonable parameters according to the atom type
#        2. avoid missing parameters (search ff96)
#        3. append missing parameters
#        """
#        p_H_type = ""
#        if (self.name, old_atom.name) not in [
#            ("CA", "CB"),
#            ("C", "N"),
#            ("CB", "CA"),
#            ("N", "C"),
#        ]:
#            raise Exception("method 1 only support truncations of CA-CB and C-N")
#        else:
#            if self.name == "CA":
#                # p_H_type = XXX
#                pass
#
#        return p_H_type
#
#    def get_protons(self):
#        """
#        get connected proton based on the connectivity map
#        ---------
#        check if connectivity is obtained. get if not.
#        """
#        return []
#
#    def get_lp_infos(self):
#        """
#        get lone pair based on the connectivity map
#        ---------
#        check if connectivity is obtained. get if not.
#        """
#        return []
#
#    def get_bond_end_atom(self):
#        """
#        get the connected atom that closer to CA in the *network*
#        -------
#        check if connectivity is obtained. get if not.
#        """
#        return None
#
#    def set_byDihedral(self, A2, A3, A4, value):
#        """
#        change the atom (self) coordinate by dihedral.
#        """
#        pass
#
#    def set_byAngle(self, A2, A3, value):
#        """
#        change the atom (self) coordinate by angle.
#        """
#        pass
#
#    def set_byBond(self, A2, value):
#        """
#        change the atom (self) coordinate by distance.
#        """
#        pass
#

    def build(self,
              a_id=None,
              r_id=None,
              c_id=None,
              ff="AMBER",
              forcefield="ff14SB"):
        """
        generate an output line. End with LF
        return a line str
        use self.id and parent id if not assigned
        -------
        must use after sort!
        """
        # default
        if a_id == None:
            a_id = self.id
        if r_id == None:
            r_id = self.residue_number
        if c_id == None:
            c_id = self.chain_id

        if ff == "AMBER":
            # build an amber style line here
            l_type = "{:<6}".format("ATOM")
            a_index = "{:>5d}".format(a_id)

            if forcefield == "ff14SB":
                # ff14SB by default
                # atom name TODO when deal with 2-letter element in artifical residue and ligand
                if len(self.atom_name) > 3:
                    a_name = "{:<4}".format(self.atom_name)
                else:
                    a_name = "{:<3}".format(self.atom_name)
                    a_name = " " + a_name
                r_name = "{:>3}".format(self.residue_name)
            else:
                raise Exception("Only support ff14SB atom/resiude name now")

            c_index = c_id
            r_index = "{:>4d}".format(r_id)
            x = "{:>8.3f}".format(self.x_coord)
            y = "{:>8.3f}".format(self.y_coord)
            z = "{:>8.3f}".format(self.z_coord)

        return f"{l_type}{a_index} {a_name} {r_name} {c_index}{r_index}    {x}{y}{z}  1.00  0.00"
#

    def build_oniom(self,
                    layer,
                    chrg=None,
                    cnt_info: list = None,
                    if_lig=0,
                    if_sol=0):
        """
        build line for oniom. Use element name for ligand atoms since they are mostly in QM regions.
        Gaussian use *ff96* which the atom type is corresponding to all_amino94.lib and ion94.lib in Amber distribution
        ---------
        chrg    : charge of current atom
        layer   : layer where atom in
        cnt_info: (low atom only) when current atom was a boundary atom. Provide [cnt_atom_ele, cnt_atom_label, cnt_atom_id]
        """
        cnt_flag = ""
        self.get_ele()
        if layer not in ["h", "l"]:
            raise Exception('build_oniom: please use: "h" or "l" for layer')
        if layer == "h":
            fz_flag = "0"
            ly_flag = "H"
        if layer == "l":
            fz_flag = "-1"
            ly_flag = "L"
            if cnt_info != None:
                cnt_flag = (" " + cnt_info[0] + "-" + cnt_info[1] + " " +
                            str(cnt_info[2]))

        # label: deal with N/C terminal and ligand
        if if_lig:
            G16_label = self.ele
        if if_sol:
            G16_label = G16_label_map[self.parent.name][self.name]
        if not if_sol and not if_lig:
            r = self.resi
            r1 = self.resi.chain.residues[0]
            rm1 = self.resi.chain.residues[-1]
            if r == r1 and self.name in ["H1", "H2", "H3"]:
                G16_label = "H"
            else:
                if r == rm1 and self.name == "OXT":
                    G16_label = "O2"
                else:
                    G16_label = G16_label_map[self.parent.name][self.name]

        # chrg
        if chrg == None:
            try:
                chrg = self.charge  # from prmtop
            except NameError:
                raise Exception(
                    "You need to at least provide a charge or use get_atom_charge to get one from prmtop file."
                )

        atom_label = "{:<16}".format(" " + self.ele + "-" + G16_label + "-" +
                                     str(round(chrg, 6)))
        fz_flag = "{:>2}".format(fz_flag)
        x = "{:<14.8f}".format(self.coord[0])
        y = "{:<14.8f}".format(self.coord[1])
        z = "{:<14.8f}".format(self.coord[2])

        line = (atom_label + " " + fz_flag + "   " + x + " " + y + " " + z +
                " " + ly_flag + cnt_flag + line_feed)

        return line


#    def get_around(self, rad):
#        pass
#
#    def get_ele(self):
#        """
#        get self.ele from a certain map according to the ff type
#        """
#        if self.name in Resi_Ele_map[self.ff].keys():
#            self.ele = Resi_Ele_map[self.ff][self.name]
#        else:
#            self.ele = re.match("^[A-Z][a-z]?", self.name).group()
#
#    """
#    ====
#    Special Method
#    ====
#    """
#
#    def __getattr__(self, key):
#        if key == "residue" or key == "resi":
#            return self.parent
#        else:
#            Exception("bad key: getattr error")
#
#    def __int__(self):
#        return self.id
#
#
