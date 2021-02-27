from helper import line_feed

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
    

    def build(self, ff='AMBER'):
        '''
        build a pdb line str
        '''
        self.get_alternate_location_indicator()
        self.get_charge()
        self.get_element()
        self.get_insert_code()
        self.get_occupancy()
        self.get_seg_id()
        self.get_temp_factor()

        if ff == 'AMBER':
            #build an amber style line here
            l_type = '{:<6}'.format(self.line_type)
            a_index = '{:>5d}'.format(self.atom_id)

            if len(self.atom_name) > 3:
                a_name = '{:<4}'.format(self.atom_name)
            else:
                a_name = '{:<3}'.format(self.atom_name)
                a_name = ' '+a_name
            r_name = '{:>3}'.format(self.resi_name)

            c_index = self.chain_id
            r_index = '{:>4d}'.format(self.resi_id)
            x = '{:>8.3f}'.format(self.atom_x)
            y = '{:>8.3f}'.format(self.atom_y)
            z = '{:>8.3f}'.format(self.atom_z)
            
            AL_id = '{:1}'.format(self.AL_id)
            insert_code = '{:1}'.format(self.insert_code)
            occupancy = '{:>6.2f}'.format(self.occupancy)
            temp_factor = '{:>6.2f}'.format(self.temp_factor)
            seg_id = '{:<4}'.format(self.seg_id)
            element = '{:1}'.format(self.element)
            charge = '{:2}'.format(self.charge)

        #example: ATOM   5350  HB2 PRO   347      32.611  15.301  24.034  1.00  0.00
        line = l_type + a_index +' '+a_name + AL_id + r_name+' '+ c_index + r_index + insert_code + '   ' + x + y + z + occupancy + temp_factor + '      ' + seg_id + element + charge + line_feed
        return line


    #misc
    def get_alternate_location_indicator(self):
        self.AL_id = self.line[16:17].strip()
        return self.AL_id
    def get_insert_code(self):
        self.insert_code = self.line[26:27].strip()
        return self.insert_code
    def get_occupancy(self):
        self.occupancy = self.line[54:60].strip()
        if self.occupancy == '':
            self.occupancy = 1.0
        else:
            self.occupancy = float(self.occupancy)
        return self.occupancy
    def get_temp_factor(self):
        self.temp_factor = self.line[60:66].strip()
        if self.temp_factor == '':
            self.temp_factor = 0.0
        else:
            self.temp_factor = float(self.temp_factor)
        return self.temp_factor
    def get_seg_id(self):
        self.seg_id = self.line[72:76].strip()
        return self.seg_id
    def get_element(self):
        self.element = self.line[76:78].strip()
        return self.element
    def get_charge(self):
        self.charge = self.line[78:80].strip()
        return self.charge
   

        
