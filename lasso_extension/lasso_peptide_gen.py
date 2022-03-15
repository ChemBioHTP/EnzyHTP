""" Extends lasso peptide from input coordinate file and writes to outfile.
"""

import numpy as np
from openbabel import pybel


def read_scaffold(filename: str) -> list:
    """Reads pdb file containing scaffold structure.

    Args:
        filename: str containing name/path of file to be read

    Return:
        Coordinate of N-Methyl attachment point from first line of scaffold file
        list containing [element_symbol, np.array(x, y, z)] for each atom

    Raises:
        ValueError: A file format other than pdb is used
    """

    if not filename.endswith(".pdb"):
        raise ValueError("Filename must end in .xyz or .pdb")

    coords = []

    with open(filename, encoding="utf-8") as infile:
        flag = int(infile.readline().split("=")[1])
        lines = infile.readlines()
        for line in lines:
            line_split = line.split()
            if len(line_split) > 1:
                symbol = line_split[2][0]
                xyz = np.array(list(map(float, line_split[5:8])))
                coords.append([symbol, xyz])

    return flag, coords


def read_extender(filename: str) -> list:
    """Reads xyz file containing extender structure.

    Args:
        filename: str containing name/path of file to be read

    Return:
        list containing [element_symbol, np.array(x, y, z)] for each atom

    Raises:
        ValueError: A file format other than xyz is used
    """

    if not filename.endswith(".xyz"):
        raise ValueError("Filename must end in .xyz")

    coords = []

    with open(filename, encoding="utf-8") as infile:
        lines = infile.readlines()[2:]
        for line in lines:
            line_split = line.split()
            symbol = line_split[0]
            xyz = np.array(list(map(float, line_split[1:])))
            coords.append([symbol, xyz])

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


def gen_xyz(coords: list, outfile_name: str) -> str:
    """Writes the calculated coordinates to an xyz-formatted string.

    Args:
        coords: The list of coordinates in [[element_symbol, xyz_coord_array], ...] form
        outfile_name: The title entry for the xyz string

    Returns:
        A string containing the coordinates in xyz format
    """

    xyz_str = str(len(coords)) + "\n"
    xyz_str += outfile_name + "\n"
    for atom in coords:
        coord_string = atom[0] + "  " if len(atom[0]) == 1 else ""
        coord_string += "".join([" " * (7 if atom[1][i] < 0 else 8) + f"{atom[1][i]:.5f}" for i in range(3)])
        xyz_str += coord_string + "\n"

    return xyz_str


def convert_to_PDB(xyz_str: str, outfile_name: str) -> None:
    """Converts an xyz-formatted string to a PDB file in the output_structures folder
    
    Args:
        xyz_str: The xyz-formatted string to read coordinates from
        outfile_name:The name of the final PDB file    
    """

    extended = pybel.readstring("xyz", xyz_str)
    extended.write("pdb", "lasso_extension/output_structures/" + outfile_name + ".pdb", overwrite=True)


def lasso_extender(scaffold_file: str,
                   tail_length: int,
                   outfile: str,
                   extension_file: str = "lasso_extension/extenders/alanine_planar_N.xyz") -> None:
    """Extends the scaffold lasso peptide using the extension amino acid by tail_length.
    
    Args:
        scaffold_file: The name of the scaffold xyz or pdb structure to add a tail to
        tail_length: The number of amino acids to add
        outfile: The name of the PDB file to write to
        extension_file: A str containing the amino acid structure to add. Defaults to alanine
        
    Returns:
        None
    """

    attachment_location, scaffold = read_scaffold(scaffold_file)
    attachment_point = scaffold[attachment_location - 1][1]

    extender = read_extender(extension_file)
    extender_vectors = [[extender[i][0], find_vector(extender, 11, i)] for i in range(12)]

    vec_extender = find_vector(extender, 12, 11)
    vec_scaffold = find_vector(scaffold, attachment_location - 3, attachment_location - 1)
    extension_matrix = rotation_matrix_from_vectors(vec_scaffold, vec_extender)

    rotated_extension = list(map(lambda vec: [vec[0], np.dot(vec[1], extension_matrix)], extender_vectors))

    scaffold = scaffold[:attachment_location - 1] + scaffold[attachment_location + 5:]
    for vector in rotated_extension:
        scaffold.append([vector[0], vector[1] + attachment_point])

    for _ in range(tail_length - 1):
        attachment_point = scaffold[-4][1]

        vec_scaffold = find_vector(scaffold, -6, -4)  # C-O single bond in alanine
        extension_matrix = rotation_matrix_from_vectors(vec_scaffold, vec_extender)

        rotated_extension = list(map(lambda vec: [vec[0], np.dot(vec[1], extension_matrix)], extender_vectors))

        scaffold = scaffold[:-4] + scaffold[-2:]  # -4 and -3 are -OH group
        for vector in rotated_extension:
            scaffold.append([vector[0], vector[1] + attachment_point])

    if outfile.endswith(".pdb"):
        outfile = outfile[:-4]
    xyz_str = gen_xyz(scaffold, outfile)

    convert_to_PDB(xyz_str, outfile)


def lasso_peptide_gen(ring: int,
                      loop: int,
                      tail: int,
                      isopeptide: str,
                      outfile: str,
                      extender: str = "lasso_extension/extenders/alanine_planar_N.xyz"):
    """Creates lasso peptide PDB with desired ring, loop, and tail size

    Args:
        ring: The number of residues in the lasso peptide ring
        loop: The number of residues in the lasso peptide upper loop
        tail: The number of residues in the lasso peptide tail
        isopeptide: The type of linkage closing the ring. asx for aspartate or glx for glutamate
        outfile: The name of the final PDB to be placed into the output_structure file
        extender: The xyz file containing the amino acid extender if not alanine
    """

    if ring != 7:
        raise ValueError("That ring size is not currently available")
    elif loop not in [7, 8]:
        raise ValueError("That loop size is not currently available")
    elif tail < 0:
        raise ValueError("Tail size must be >= 0")
    elif isopeptide not in ["asx", "glx"]:
        raise ValueError("Isopeptide linker must be asx (aspartate) or glx (glutamate)")

    scaffold = f"lasso_extension/scaffolds/{ring}mr_lassoPeptide_{isopeptide}_{loop}A.pdb"
    lasso_extender(scaffold, tail, outfile, extender)


def main():
    """Driver for lasso_peptide_gen()

    """

    ring = 7
    loop = 7
    tail_length = 4
    isopeptide = "asx"
    outfile = "asx_7A"

    lasso_peptide_gen(ring, loop, tail_length, isopeptide, outfile)


if __name__ == "__main__":
    main()
