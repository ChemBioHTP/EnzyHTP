'''
Extends lasso peptide from xyz file
'''

import numpy as np

def read_structure(filename : str) -> list:
    '''
    Reads xyz file containing structure
    Returns list containing [element_symbol, x, y, z] for each atom
    '''
    coords = []
    with open(filename, encoding="utf-8") as infile:
        lines = infile.readlines()[2:]
        for line in lines:
            line_split = line.split()
            symbol = line_split[0]
            xyz = np.array(list(map(float, line_split[1:])))
            coords.append([symbol, xyz])
    return coords

def find_vector(structure: list, coord_1 : int, coord_2 : int) -> np.array:
    '''
    Finds vector between atoms in the structure of the given coordinates
    The vector points from coord_1 to coord_2
    coord_1 and coord_2 use 0-based indexing
    '''
    return structure[coord_2][1] - structure[coord_1][1]

def rotation_matrix_from_vectors(vec1 : np.array, vec2 : np.array) -> np.matrix:
    '''
    Calculates rotation matrix between two vectors which must be non-parallel
    Code borrowed from:
    https://stackoverflow.com/questions/63525482/finding-the-rotation-matrix-between-two-vectors-in-python
    '''

    if type(vec1) != np.ndarray or type(vec2) != np.ndarray:
        raise ValueError("Input vectors must be numpy arrays")
    if vec1.shape != (3,) or vec2.shape != (3,):
        raise ValueError("Input vectors must be 3D vectors")
    vec1_norm = (vec1 / np.linalg.norm(vec1)).reshape(3)
    vec2_norm = (vec2 / np.linalg.norm(vec2)).reshape(3)
    
    cross = np.cross(vec1_norm, vec2_norm)
    dot = np.dot(vec1_norm, vec2_norm)
    cross_mag = np.linalg.norm(cross)
    if cross_mag == 0:
        return np.eye(3)
    kmat = np.array([[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - dot) / (cross_mag ** 2))
    return rotation_matrix

def write_xyz(coords : list, outfile_name : str) -> None:
    '''
    Writes the calculated coordinates to an xyz file for further analysis
    '''
    if not outfile_name.endswith(".xyz"):
        outfile_name += ".xyz"
    
    with open("lasso_extension/structures/"+outfile_name, "w+", encoding="utf-8") as outfile:
        print(len(coords), file=outfile)
        print(outfile_name, file=outfile)

        for atom in coords:
            coord_string = atom[0] + "  " if len(atom[0]) == 1 else "" 
            coord_string += "".join([" " * (7 if atom[1][i] < 0 else 8) +f"{atom[1][i]:.5f}" for i in range(3)])
            print(coord_string, file=outfile)
        

def main():
    '''
    Reads in input PDB file and amino acid to extend
    Concatenates tail_len instances of that amino acid to the lasso peptide
    '''
    tail_len = 4
    
    # TO/DO: add input parameters (command line args or function parameters -- talk to team)
    extender = read_structure("lasso_extension/structures/alanine_planar_N.xyz")
    base = read_structure("lasso_extension/structures/severedLassoPeptide.xyz")
   
    extender_vectors = [[extender[i][0], find_vector(extender, 11, i)] for i in range(12)]
    
    # 1: determine vectors vec_extender, vec_base, and attachment_point
    vec_extender = find_vector(extender, 12, 11)
    vec_base = find_vector(base, 27, 35)

    attachment_point = base[35][1]

    # 2: find matrix between vec_extender and vec_base
    extension_matrix = rotation_matrix_from_vectors(vec_base, vec_extender)

    
    # 3: apply matrix to all atoms in extender
    rotated_extension = list(map(lambda vec: [vec[0], np.dot(vec[1], extension_matrix)], extender_vectors))
    
    # 4: trim hanging group from base (H3C-N-H moiety or OH)
    base = base[:35]+base[41:] # N36-H41 are H3C-N-H group

    # 5: recombine rotated extension at attachment_point
    for vector in rotated_extension:
        base.append([vector[0], vector[1] + attachment_point])
    

    # 6: go back to step 1 until desired length is reached
    for _ in range(tail_len - 1):
        vec_base = find_vector(base, -6, -4) # C-O bond in alanine

        attachment_point = base[-4][1]

        # 2: find matrix between vec_extender and vec_base
        extension_matrix = rotation_matrix_from_vectors(vec_base, vec_extender)

        
        # 3: apply matrix to all atoms in extender
        rotated_extension = list(map(lambda vec: [vec[0], np.dot(vec[1], extension_matrix)], extender_vectors))
        
        # 4: trim hanging group from base (H3C-N-H moiety or -OH)
        base = base[:-4]+base[-2:] # -4 and -3 are -OH group

        # 5: recombine rotated extension at attachment_point
        for vector in rotated_extension:
            base.append([vector[0], vector[1] + attachment_point])
        

    # 7: Write to output file
    write_xyz(base, "extendedPeptideMulti")



if __name__ == "__main__":
    main()
