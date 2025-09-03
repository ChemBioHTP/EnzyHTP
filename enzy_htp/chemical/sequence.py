"""Sequence utilities for handling protein and nucleic acid sequences.

This module provides helper functions for creating FASTA files, validating sequences,
and other sequence-related operations used throughout EnzyHTP.

Author: Qianzhen (QZ) Shao <shaoqz@icloud.com>
Date: 2025-08-18
"""
from __future__ import annotations
import tempfile
from typing import List, Union, Tuple
from pathlib import Path


def create_fasta_from_sequences(
    sequences: Union[str, List[str]], 
    sequence_ids: Union[str, List[str], None] = None,
    output_path: Union[str, Path, None] = None,
) -> str:
    """Create a FASTA file from protein sequences.
    
    Args:
        sequences: Single sequence string or list of sequence strings
        sequence_ids: Optional sequence identifiers. If None, auto-generates seq_0, seq_1, etc.
        output_path: Optional output file path. If None, creates temporary file.
        delete_on_close: If True and output_path is None, creates temporary file that auto-deletes
        
    Returns:
        str: Path to the created FASTA file
        
    Raises:
        ValueError: If sequences and sequence_ids lists have different lengths
        TypeError: If sequences is not str or list of str
        
    Example:
        >>> fasta_path = create_fasta_from_sequences(["ACDEFG", "GHIKLM"])
        >>> # Creates file with content:
        >>> # >seq_0
        >>> # ACDEFG
        >>> # >seq_1
        >>> # GHIKLM
    """
    # Normalize inputs
    if isinstance(sequences, str):
        sequences = [sequences]
    if not isinstance(sequences, (list, tuple)):
        raise TypeError("sequences must be str or list/tuple of str")
    
    if sequence_ids is None:
        sequence_ids = [f"seq_{i}" for i in range(len(sequences))]
    elif isinstance(sequence_ids, str):
        sequence_ids = [sequence_ids]
    
    if len(sequences) != len(sequence_ids):
        raise ValueError(f"sequences ({len(sequences)}) and sequence_ids ({len(sequence_ids)}) must have same length")
    
    # Create output file
    if output_path is None:
        temp_file = tempfile.NamedTemporaryFile(
            mode='w', 
            suffix='.fasta',
            delete=False,
        )
        file_handle = temp_file
        file_path = temp_file.name
    else:
        file_path = str(output_path)
        file_handle = open(file_path, 'w')
    
    try:
        # Write FASTA content
        for seq_id, sequence in zip(sequence_ids, sequences):
            file_handle.write(f">{seq_id}\n{sequence}\n")
        
        if output_path is not None:
            file_handle.close()
        else:
            file_handle.flush()
            
        return file_path
        
    finally:
        file_handle.close()


def parse_fasta_file(fasta_path: Union[str, Path]) -> List[Tuple[str, str]]:
    """Parse a FASTA file and return sequence ID and sequence pairs.
    
    Args:
        fasta_path: Path to the FASTA file
        
    Returns:
        List of (sequence_id, sequence) tuples
        
    Raises:
        FileNotFoundError: If FASTA file doesn't exist
        ValueError: If FASTA file is malformed
        
    Example:
        >>> sequences = parse_fasta_file("sequences.fasta")
        >>> # Returns: [("seq_0", "ACDEFG"), ("seq_1", "GHIKLM")]
    """
    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    sequences = []
    current_id = None
    current_seq = []
    
    with open(fasta_path, 'r') as file:
        for line_num, line in enumerate(file, 1):
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_id is not None:
                    sequences.append((current_id, ''.join(current_seq)))
                
                # Start new sequence
                current_id = line[1:].strip()
                current_seq = []
                
            elif current_id is not None:
                current_seq.append(line)
                
            else:
                raise ValueError(f"Malformed FASTA file at line {line_num}: sequence data before header")
    
    # Save last sequence
    if current_id is not None:
        sequences.append((current_id, ''.join(current_seq)))
    
    if not sequences:
        raise ValueError("No sequences found in FASTA file")
    
    return sequences


def validate_protein_sequence(sequence: str, allow_ambiguous: bool = True) -> bool:
    """Validate if a string is a valid protein sequence.
    
    Args:
        sequence: Protein sequence string to validate
        allow_ambiguous: If True, allows ambiguous amino acids (B, J, O, U, X, Z)
        
    Returns:
        bool: True if sequence is valid, False otherwise
        
    Example:
        >>> validate_protein_sequence("ACDEFGHIKLMNPQRSTVWY")
        True
        >>> validate_protein_sequence("ACDEFGHIKLMNPQRSTVWYXZ")
        True
        >>> validate_protein_sequence("ACDEFGHIKLMNPQRSTVWYXZ", allow_ambiguous=False)
        False
    """
    if not isinstance(sequence, str):
        return False
    
    # Standard 20 amino acids
    standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
    
    # Ambiguous amino acids
    ambiguous_aa = set("BJOUX Z")  # B=D/N, J=I/L, O=Pyrrolysine, U=Selenocysteine, X=any, Z=E/Q
    
    valid_chars = standard_aa
    if allow_ambiguous:
        valid_chars = valid_chars.union(ambiguous_aa)
    
    sequence_chars = set(sequence.upper())
    return sequence_chars.issubset(valid_chars)


def get_sequence_length(sequence: str) -> int:
    """Get the length of a protein sequence, ignoring whitespace and newlines.
    
    Args:
        sequence: Protein sequence string
        
    Returns:
        int: Length of the cleaned sequence
        
    Example:
        >>> get_sequence_length("ACDE FGH\\nIKL")
        11
    """
    return len(''.join(sequence.split()))


def clean_sequence(sequence: str) -> str:
    """Clean a protein sequence by removing whitespace and converting to uppercase.
    
    Args:
        sequence: Raw protein sequence string
        
    Returns:
        str: Cleaned sequence string
        
    Example:
        >>> clean_sequence("acde fgh\\nikl")
        "ACDEFGHIKL"
    """
    return ''.join(sequence.split()).upper()