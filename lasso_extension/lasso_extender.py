""" Extends lasso peptide from input coordinate file and writes to outfile.
"""

import numpy as np


def read_structure(filename: str) -> list:
    """Reads xyz or pdb file containing structure.

    Args:
        filename: str containing name/path of file to be read

    Return:
        list containing [element_symbol, np.array(x, y, z)] for each atom

    Raises:
        ValueError: A file format other than xyz or pdb is used
    """
    coords = []

    with open(filename, encoding="utf-8") as infile:
        if filename.endswith(".xyz"):
            lines = infile.readlines()[2:]
            for line in lines:
                line_split = line.split()
                symbol = line_split[0]
                xyz = np.array(list(map(float, line_split[1:])))
                coords.append([symbol, xyz])
        elif filename.endswith(".pdb"):
            lines = infile.readlines()
            for line in lines:
                line_split = line.split()
                if len(line_split) > 1:
                    symbol = line_split[-1]
                    xyz = np.array(list(map(float, line_split[5:8])))
                    coords.append([symbol, xyz])
        else:
            infile.close()
            raise ValueError("Filename must end in .xyz or .pdb")
    return coords


def find_vector(structure: list, coord_1: int, coord_2: int) -> np.array:
    """Finds vector between atoms in the structure of the given coordinates.
    
    Args:
        structure: The list of [[element_symbol, coordinate_array],...]
        coord_1: The index (0-based) of the origin of the vector
        coord_2: The index (0-based) of the tip of the vector
    
    Returns:
        The vector points from coord_1 to coord_2
    """
    return structure[coord_2][1] - structure[coord_1][1]


def rotation_matrix_from_vectors(vec1: np.array, vec2: np.array) -> np.matrix:
    """Calculates rotation matrix between two vectors.

    Code borrowed from:
    https://stackoverflow.com/questions/63525482/finding-the-rotation-matrix-between-two-vectors-in-python

    Args:
        vec1: The base vector
        vec2: The target vector
    
    Returns:
        A 3 by 3 numpy matrix that represents the rotational transformation from vec 1 to vec 2.

    Raises:
        ValueError: if vec1 and vec2 are not numpy arrays of the proper dimension.
    """

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
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - dot) / (cross_mag**2))
    return rotation_matrix


def write_xyz(coords: list, outfile_name: str) -> None:
    """Writes the calculated coordinates to an xyz file.

    Args:
        coords: The list of coordinates in [[element_symbol, xyz_coord_array], ...] form
        outfile_name: The name of the xyz file to write the coordinates

    Returns:
        None
    """
    if not outfile_name.endswith(".xyz"):
        outfile_name += ".xyz"

    with open("lasso_extension/structures/" + outfile_name, "w+", encoding="utf-8") as outfile:
        print(len(coords), file=outfile)
        print(outfile_name, file=outfile)

        for atom in coords:
            coord_string = atom[0] + "  " if len(atom[0]) == 1 else ""
            coord_string += "".join([" " * (7 if atom[1][i] < 0 else 8) + f"{atom[1][i]:.5f}" for i in range(3)])
            print(coord_string, file=outfile)


def lasso_extender(base_file: str,
                   tail_length: int,
                   outfile: str,
                   extension_file: str = "lasso_extension/structures/alanine_planar_N.xyz") -> None:
    """Extends the base lasso peptide using the extension amino acid by tail_length.
    
    Args:
        base_file: The name of the base xyz or pdb structure to add a tail to
        tail_length: The number of amino acids to add
        outfile: The name of the file to write to 
        extension_file: A str containing the amino acid structure to add. Defaults to alanine
        
    Returns:
        None
    """

    base = read_structure(base_file)
    attachment_point = base[35][1]

    extender = read_structure(extension_file)
    extender_vectors = [[extender[i][0], find_vector(extender, 11, i)] for i in range(12)]

    vec_extender = find_vector(extender, 12, 11)
    vec_base = find_vector(base, 27, 35)
    extension_matrix = rotation_matrix_from_vectors(vec_base, vec_extender)

    rotated_extension = list(map(lambda vec: [vec[0], np.dot(vec[1], extension_matrix)], extender_vectors))

    base = base[:35] + base[41:]  # N36-H41 are H3C-N-H group
    for vector in rotated_extension:
        base.append([vector[0], vector[1] + attachment_point])

    for _ in range(tail_length - 1):
        attachment_point = base[-4][1]

        vec_base = find_vector(base, -6, -4)  # C-O single bond in alanine
        extension_matrix = rotation_matrix_from_vectors(vec_base, vec_extender)

        rotated_extension = list(map(lambda vec: [vec[0], np.dot(vec[1], extension_matrix)], extender_vectors))

        base = base[:-4] + base[-2:]  # -4 and -3 are -OH group
        for vector in rotated_extension:
            base.append([vector[0], vector[1] + attachment_point])

    write_xyz(base, outfile)


def main():
    """Driver for lasso_extender()

    Extends the severedLassoPeptide.xyz file by 4 of the default Alanine extender
    """
    base = "lasso_extension/structures/severedLassoPeptide.xyz"
    tail_length = 4
    outfile = "extendedPeptideMulti"
    lasso_extender(base, tail_length, outfile)


if __name__ == "__main__":
    main()
