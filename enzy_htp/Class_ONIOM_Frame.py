'''
- change/update coordinate in a ONIOM template file
- need all files in operation share a **same** atomic/coordinate order !!
Usage:
1. Obtain the coordinate
    - from .mdcrd file:  coords = Frame.fromMDCrd(path)  // return a list of all frames in the mdcrd file. You may want to use cpptraj to sample the wanted frame into the file
    - from Gaussian output file: coord = Frame.fromGaussinOut(path) // return the last point of the gaussian opt/freq
2. (optional) shift some orders of the coordinate
    - frame.shift_line(shift_list) // shift_list is a list of (l1, l2): l1 is the moving line, l2 is the line before the target position. *l2 cannot be same as any l1 in the list.
3. combine the coordinate with the template
    - write_to_template(template_path, frame_obj) // Use out_path and index to customize the output filename and path.
-----------------
4. select and write a visible file
    select:
    - Frame.sele_unfreeze(g_file)
    - Frame.sele_high(g_file) // return a dict of select infomation (index : atom_name)
    write:
    - frame.write_sele_lines(self, sele_list, out_path='sele_coord.'+ff, ff='gjf')
-----------------
5. Read Frequencies
    - getFreq(g_out_file) // return a list of frequence numbers
-----------------
'''
import numpy as np
from Class_Conf import Config
from helper import line_feed, set_distance
import re
import os

# In gjf: 
#   pattern for determining the beginning of the coordinate (strip)
ChrgSpin_pattern = r'(?: ?\-?\+?[0-9] [0-9])+'
#   pattern for unfreezed atoms
unfreeze_pattern = r'[A-z,\-,0-9,\.]+ +0 '
#   pattern for high layer atoms
high_pattern = r'[0-9]+ +H'
# In mdcrd: 
#   pattern for determining the end of a old frame (no strip)
digit_pattern = r'[ ,\-,0-9][ ,\-,0-9][ ,\-,0-9][0-9]\.[0-9][0-9][0-9]'
frame_sep_pattern = digit_pattern * 3 + line_feed
# In log/out:
#   pattern for determining the position of frequencies
freq_pattern = r'Frequencies'





class Frame:
    '''
    store and operate massive coordinate for ONIOM I/O
    ---------
    raw coordniate only. use .prmtop to relate to the chemistry info 
    '''

    def __init__(self, coord):
        '''
        coord: a 2D list of coordinate of each atom [[x,y,z],...]
        '''
        self.coord = coord

    @classmethod
    def fromMDCrd(cls, mdcrd_file):
        '''
        read a list of coordinates for frames
        -------
        return a list of Frame object
        '''
        coords = []
        coord = []
        atom_coord = []
        counter = 1
        end_flag_1 = 0
        fake_end_flag = 0
        #os.system('echo >> '+mdcrd_file)
        debug_counter = 0

        with open(mdcrd_file) as f:
            while True:
                # use the while True format to detect the EOF
                line=f.readline()            

                if re.match(digit_pattern, line) == None:
                    # not data line // could be faster here
                    if not end_flag_1 and not fake_end_flag:
                        continue
                
                if fake_end_flag:
                    # if next line of fake end is the file end 
                    if line == '':
                        break

                if end_flag_1:
                    debug_counter += 1
                    if re.match(frame_sep_pattern, line) != None:
                        # last line is a fake end line
                        coord.append(holder)
                        coords.append(cls(coord))
                        # empty for next loop
                        coord = []
                        end_flag_1 = 0
                        fake_end_flag = 1
                        continue
                    else:                        
                        # last line is a real end line
                        coords.append(cls(coord))
                        # empty for next loop
                        coord = []
                        end_flag_1 = 0
                        if line == '':
                            #the last line
                            break
                        if line == line_feed:
                            if Config.debug >= 1:
                                print("Frame.fromMDCrd: WARNING: unexpected empty line detected. Treat as EOF. exit reading")
                            break
                        # do not skip if normal next line

                else:
                    if re.match(frame_sep_pattern, line) != None:
                        # find line that mark the end of a frame // possible fake line
                        # store the last frame and empty the holder
                        end_flag_1 = 1 
                        lp = re.split(' +', line.strip())
                        # hold the info
                        holder = [float(lp[0]), float(lp[1]), float(lp[2])]
                        continue
                

                # normal data lines
                lp = re.split(' +', line.strip())
                for i in lp:
                    if counter < 3:
                        atom_coord.append(float(i))
                        counter = counter + 1
                    else:
                        atom_coord.append(float(i))
                        coord.append(atom_coord)
                        # empty for next atom
                        atom_coord = []        
                        counter = 1
        return coords

    @classmethod
    def fromGaussinOut(cls, g_out_file):
        '''
        get last step from the Gaussian out file, according to the Input orientation
        '''
        coord = []
        f_counter = 0
        f_counter_2 = 0
        coord_flag = 0
        skip4 = 0
        with open(g_out_file) as f:
            # get where is the last 'Input orientation:' mark by count the number
            for line in f:
                if line.strip() == 'Input orientation:':
                    f_counter = f_counter + 1

        # reopen to avoid the mem problem
        with open(g_out_file) as f:
            # get the last part based on the count
            for line in f:
                if line.strip() == 'Input orientation:':
                    f_counter_2 = f_counter_2 + 1
                    if f_counter_2 == f_counter:
                        #here is the last frame
                        coord_flag = 1
                
                if coord_flag:
                    # skip 4 title lines
                    if skip4 <= 4:
                        skip4 = skip4 + 1
                        continue
                    # detect end of the coord section
                    if line.strip() == '---------------------------------------------------------------------':
                        coord_flag = 0
                        break

                    lp = line.strip().split()
                    atom_coord = [float(lp[3]), float(lp[4]), float(lp[5])]
                    coord.append(atom_coord)
        
        return cls(coord)


    def shift_line(self, shift_list):
        '''
        shift_list = [(l1,l2), ... ]
        -------------
        shift l1 to the place **after** l2
        l2 must not be any other l1 in the list
        follow the same order in the shift_list when there's more than one same l2. (e.g. [(5,1), (6,1)])
        '''
        
        # convert to a map (memorize the old line order) & capture target lines: (l2, l1_coord)
        new_coord = {}
        t_lines = []
        for i, atom_coord in enumerate(self.coord):
            if_t_line = 0
            for l1,l2 in shift_list:
                if i+1 == l1:
                    t_lines.append((str(l2),atom_coord))
                    if_t_line = 1
            if not if_t_line:        
                new_coord[str(i+1)] = atom_coord
        
        # reorder and save the new coordinate
        self.coord = []
        for key in new_coord.keys():
            self.coord.append(new_coord[key])
            for l2, l1_coord in t_lines:
                if key == l2:                    
                    self.coord.append(l1_coord)


    def write_to_template(self, t_file_path, out_path=None, index : str = None, ifchk=1):
        '''
        1. find the beginning and ending of the coordinate section.
        2. replace coordinate based on the same atom sequence.
        ------
        write to out_path
        '''
        if out_path == None:
            out_path = t_file_path[:-4]+'_newcoord.gjf'
        if index != None:
            out_path = out_path[:-4]+'_'+index+'.gjf'
        chk_path = out_path[:-3]+'chk'
        
        with open(t_file_path) as f:
            with open(out_path,'w') as of:
                coord_b_flag = 0
                coord_e_flag = 0
                coord_index = 0
                for line in f:
                    
                    if re.search('chk_place_holder', line.strip()) != None:
                        if ifchk:
                            of.write(r'%chk='+chk_path + line_feed)
                            continue
                        else:
                            continue

                    if re.match(ChrgSpin_pattern, line.strip()) != None:
                        coord_b_flag = 1
                        of.write(line)
                        continue
                    
                    # end when a space line appears
                    if coord_b_flag:
                        if line == line_feed:
                            coord_e_flag = 1

                    if coord_b_flag and not coord_e_flag:
                        # line part
                        lp = line.strip().split()
                        line_coord = self.coord[coord_index]
                        
                        # potential *CHANGE* here about where coord is in a line
                        label = '{:<20}'.format(lp[0])
                        freeze_mark = '{:<3}'.format(lp[1])
                        x = '{:>15.8f}'.format(line_coord[0])
                        y = '{:>15.8f}'.format(line_coord[1])
                        z = '{:>15.8f}'.format(line_coord[2])
                        layer_mark = lp[5]

                        new_line = ' '+label+' '+freeze_mark+x+y+z+' '+layer_mark

                        if len(lp) > 6:
                            new_line = new_line + ' ' + lp[6] + ' ' + lp[7]
                        new_line = new_line + line_feed

                        of.write(new_line)
                        coord_index = coord_index+1
                        continue
                    
                    of.write(line)
                    
        if ifchk:
            return out_path, chk_path

        return out_path


    def write_sele_lines(self, sele_list, out_path = None, ff='gjf', g_route=None, chrgspin=None, ifchk=0):
        '''
        the sele list should be a map like: {'int':'atom_name', '1':'H', ....}
        --------
        obtain atomic info externally. Mimic the Amber prmtop strategy ---> only have the minimium atomic info
        TODO: Any meaning to require a map form of sele_list? instead of a list.
        '''
        if out_path == None:
            out_path ='sele_coord.'+ff
        if ifchk:
            chk_path = out_path[:-len(ff)]+'chk'

        sele_lines = []
        for sele in sele_list.keys():
            # decode val fix atoms
            if '-' in sele:
                ele_name = sele_list[sele]
                fix_info = sele.split('-')
                fix_val_coord = self._get_fix_val_coord(fix_info[0], fix_info[1], fix_info[2])
                sele_lines.append((ele_name, fix_val_coord))
                continue

            # get coord for normal atoms
            # clean up
            if sele[-1] not in '1234567890':
                sele_id = sele[:-1]
            else:
                sele_id = sele 
            # get coord
            sele_lines.append((sele_list[sele], self.coord[int(sele_id)-1]))
        
        with open(out_path, 'w') as of:
            if ff != 'xyz' and ff != 'gjf':
                raise Exception('only support xyz and gjf format now')
            if ff == 'xyz':
                of.write(str(len(sele_list.keys()))+line_feed)
            if ff == 'gjf':
                if ifchk:
                    of.write(r'%chk='+chk_path+line_feed)
                of.write('%mem='+str(Config.max_core*Config.n_cores)+'MB'+line_feed)
                of.write('%nprocshared='+str(Config.n_cores)+line_feed)
                if g_route == None:
                    of.write('# hf/3-21g'+line_feed)
                else:
                    of.write(g_route+line_feed)
                of.write(line_feed)
                of.write('Title Card Required'+line_feed)
                of.write(line_feed)
                if chrgspin == None:
                    of.write('0 1'+line_feed)
                else:
                    chrg = str(chrgspin[0])
                    spin = str(chrgspin[1])
                    of.write(chrg+' '+spin+line_feed)

            for atom, line_coord in sele_lines:
                label = '{:<5}'.format(atom)
                x = '{:>15.8f}'.format(line_coord[0])
                y = '{:>15.8f}'.format(line_coord[1])
                z = '{:>15.8f}'.format(line_coord[2])

                of.write(label+x+y+z+line_feed)
            # write a blank line to support "g16 < .gjf > .out" mode of gaussian
            of.write(line_feed)

    @classmethod
    def sele_unfreeze(cls, g_file, ff='gjf'):
        '''
        select unfreeze part from the gaussian oniom I/O file.
        ---------
        recommand using oniom gjf file
        '''
        sele_list = {}
        coord_index = None
        if ff == 'gjf':
            with open(g_file) as f:
                for line in f:
                    if re.match(ChrgSpin_pattern, line) != None:
                        coord_index = 0
                        continue
                    if coord_index != None:
                        coord_index = coord_index + 1 #current index

                    if re.match(unfreeze_pattern, line.strip()) != None:
                        lp = line.strip().split()
                        atom = lp[0].split('-')[0]
                        sele_list[str(coord_index)] = atom
            return sele_list

        raise Exception('only support gjf now')

    @classmethod
    def sele_high(cls, g_file, ff='gjf'):
        '''
        select high layer from the gaussian oniom I/O file.
        ---------
        recommand using oniom gjf file
        '''
        sele_list = {}
        coord_index = None
        if ff == 'gjf':
            with open(g_file) as f:
                for line in f:
                    if re.match(ChrgSpin_pattern, line) != None:
                        coord_index = 0
                        continue
                    if coord_index != None:
                        coord_index = coord_index + 1 #current index

                    if re.search(high_pattern, line.strip()) != None:
                        lp = line.strip().split()
                        atom = lp[0].split('-')[0]
                        sele_list[str(coord_index)] = atom
            return sele_list

        raise Exception('only support gjf now')
                

    def _get_fix_val_coord(self, id1, id2, d):
        '''
        search for coord of id1 and id2 and generate a new coord using d for center atom of the fix
        '''
        # search for p1
        for index, line in enumerate(self.coord):            
            if int(id1) == index + 1:
                p1 = line
                break
        for index, line in enumerate(self.coord):            
            if int(id2) == index + 1:
                p2 = line
                break
        d = float(d)
        fix_val_coord = set_distance(p1,p2,d)

        return fix_val_coord
                
    '''
    special method
    '''
    def __getitem__(self, key: int):
        '''
        Frame_obj[int]: frame.coord[int]
        '''
        return self.coord[key]



def getFreq(g_out_file):
    '''
    Get frequencies from a gaussian output file.
    --------
    return a list of frequencies. (value) 
    '''
    freqs = []
    with open(g_out_file) as f:
        for line in f:
            if re.match(freq_pattern, line.strip()):
                lp = line.strip().split()
                lp = lp[2:]

                for freq in lp:
                    freqs.append(float(freq.strip()))
    
    return freqs

