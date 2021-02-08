class PDB_line(object):
    '''
    Class for decoding the line in the PDB file. Functions used internally. 
    Generate from:
    line: PDB_line(line) --> PDB_line
    lines: PDB_line.fromlines(lines) --> [PDB_line, ...]
    '''


    def __init__(self,line):
        '''
        initilize with a specific line in the PDB file.
        Get self.line_type
        '''
        try:
            self.line=line
            self.line_type = self.line[0:6].strip()
            # atom
            self.atom_id = int(self.line[6:11])
            self.atom_name = self.line[12:16].strip()
            # residue
            self.resi_name = self.line[17:20].strip()
            self.resi_id = int(self.line[22:26])
            # chain
            self.chain_id = self.line[21:22]
            # coord
            self.atom_x = float(self.line[30:38])
            self.atom_y = float(self.line[38:46])
            self.atom_z = float(self.line[46:54])
        except ValueError:
            pass
    
    @classmethod
    def fromlines(cls, lines):
        '''
        capable with multiple lines
        return a list of PDB_line object
        '''
        line_feed='\n' # potential update

        line_list = lines.strip().split(line_feed)
        pdb_l_list = []
        for line in line_list:
            pdb_l_list.append(cls(line))

        return pdb_l_list

    #misc
    def get_alternate_location_indicator(self):
        self.AL_id = self.line[16:17].strip()
    def get_insert_code(self):
        self.insert_code = self.line[26:27]
    def get_occupancy(self):
        self.occupancy = self.line[54:60].strip()
    def get_temp_factor(self):
        self.temp_factor = self.line[60:66].strip()
    def get_seg_id(self):
        self.seg_id = self.line[72:76].strip()
    def get_element(self):
        self.element = self.line[76:78].strip()
   

        
