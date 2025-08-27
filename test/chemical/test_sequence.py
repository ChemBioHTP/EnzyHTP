"""Test module for enzy_htp.chemical.sequence

Author: Claude Code
Date: 2025-08-18
"""
import pytest
import tempfile
import os
from pathlib import Path

from enzy_htp.chemical.sequence import (
    create_fasta_from_sequences,
    parse_fasta_file,
    validate_protein_sequence,
    get_sequence_length,
    clean_sequence
)


class TestCreateFastaFromSequences:
    """Test class for create_fasta_from_sequences function."""

    def test_single_sequence_string(self):
        """Test creating FASTA from single sequence string."""
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        fasta_path = create_fasta_from_sequences(sequence)
        
        try:
            with open(fasta_path, 'r') as f:
                content = f.read()
            
            expected = ">seq_0\nACDEFGHIKLMNPQRSTVWY\n"
            assert content == expected
        finally:
            os.unlink(fasta_path)

    def test_multiple_sequences_list(self):
        """Test creating FASTA from multiple sequences."""
        sequences = ["ACDEFGHIKLMNPQRSTVWY", "DEFGHIKLMNPQRSTVWY"]
        fasta_path = create_fasta_from_sequences(sequences)
        
        try:
            with open(fasta_path, 'r') as f:
                content = f.read()
            
            expected = ">seq_0\nACDEFGHIKLMNPQRSTVWY\n>seq_1\nDEFGHIKLMNPQRSTVWY\n"
            assert content == expected
        finally:
            os.unlink(fasta_path)

    def test_custom_sequence_ids(self):
        """Test creating FASTA with custom sequence IDs."""
        sequences = ["ACDEFG", "GHIKLM"]
        seq_ids = ["protein1", "protein2"]
        fasta_path = create_fasta_from_sequences(sequences, seq_ids)
        
        try:
            with open(fasta_path, 'r') as f:
                content = f.read()
            
            expected = ">protein1\nACDEFG\n>protein2\nGHIKLM\n"
            assert content == expected
        finally:
            os.unlink(fasta_path)

    def test_custom_output_path(self):
        """Test creating FASTA with custom output path."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "test.fasta"
            sequences = ["ACDEFG"]
            
            result_path = create_fasta_from_sequences(sequences, output_path=output_path)
            
            assert result_path == str(output_path)
            assert output_path.exists()
            
            with open(output_path, 'r') as f:
                content = f.read()
            assert content == ">seq_0\nACDEFG\n"

    def test_mismatched_lengths_error(self):
        """Test error when sequences and IDs have different lengths."""
        sequences = ["ACDEFG", "GHIKLM"]
        seq_ids = ["protein1"]
        
        with pytest.raises(ValueError, match="must have same length"):
            create_fasta_from_sequences(sequences, seq_ids)

    def test_invalid_sequences_type_error(self):
        """Test error when sequences is not str or list."""
        with pytest.raises(TypeError, match="must be str or list of str"):
            create_fasta_from_sequences(123)


class TestParseFastaFile:
    """Test class for parse_fasta_file function."""

    def test_parse_simple_fasta(self):
        """Test parsing a simple FASTA file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\nACDEFG\n>seq2\nGHIKLM\n")
            fasta_path = f.name
        
        try:
            sequences = parse_fasta_file(fasta_path)
            expected = [("seq1", "ACDEFG"), ("seq2", "GHIKLM")]
            assert sequences == expected
        finally:
            os.unlink(fasta_path)

    def test_parse_multiline_sequences(self):
        """Test parsing FASTA with sequences split across multiple lines."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\nACDE\nFGHI\n>seq2\nGHIK\nLMNP\n")
            fasta_path = f.name
        
        try:
            sequences = parse_fasta_file(fasta_path)
            expected = [("seq1", "ACDEFGHI"), ("seq2", "GHIKLMNP")]
            assert sequences == expected
        finally:
            os.unlink(fasta_path)

    def test_parse_with_empty_lines(self):
        """Test parsing FASTA with empty lines."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n\nACDEFG\n\n>seq2\n\nGHIKLM\n")
            fasta_path = f.name
        
        try:
            sequences = parse_fasta_file(fasta_path)
            expected = [("seq1", "ACDEFG"), ("seq2", "GHIKLM")]
            assert sequences == expected
        finally:
            os.unlink(fasta_path)

    def test_file_not_found_error(self):
        """Test error when FASTA file doesn't exist."""
        with pytest.raises(FileNotFoundError):
            parse_fasta_file("nonexistent.fasta")

    def test_malformed_fasta_error(self):
        """Test error when FASTA file is malformed."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write("ACDEFG\n>seq1\nGHIKLM\n")  # Sequence before header
            fasta_path = f.name
        
        try:
            with pytest.raises(ValueError, match="sequence data before header"):
                parse_fasta_file(fasta_path)
        finally:
            os.unlink(fasta_path)

    def test_empty_fasta_error(self):
        """Test error when FASTA file is empty."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write("")
            fasta_path = f.name
        
        try:
            with pytest.raises(ValueError, match="No sequences found"):
                parse_fasta_file(fasta_path)
        finally:
            os.unlink(fasta_path)


class TestValidateProteinSequence:
    """Test class for validate_protein_sequence function."""

    def test_valid_standard_sequence(self):
        """Test validation of standard amino acid sequence."""
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        assert validate_protein_sequence(sequence) is True

    def test_valid_with_ambiguous_allowed(self):
        """Test validation with ambiguous amino acids allowed."""
        sequence = "ACDEFGHIKLMNPQRSTVWYXZ"
        assert validate_protein_sequence(sequence, allow_ambiguous=True) is True

    def test_invalid_with_ambiguous_disallowed(self):
        """Test validation with ambiguous amino acids disallowed."""
        sequence = "ACDEFGHIKLMNPQRSTVWYXZ"
        assert validate_protein_sequence(sequence, allow_ambiguous=False) is False

    def test_lowercase_sequence(self):
        """Test validation of lowercase sequence."""
        sequence = "acdefghiklmnpqrstvwy"
        assert validate_protein_sequence(sequence) is True

    def test_invalid_characters(self):
        """Test validation with invalid characters."""
        sequence = "ACDEFG123"
        assert validate_protein_sequence(sequence) is False

    def test_non_string_input(self):
        """Test validation with non-string input."""
        assert validate_protein_sequence(123) is False
        assert validate_protein_sequence(None) is False


class TestGetSequenceLength:
    """Test class for get_sequence_length function."""

    def test_simple_sequence(self):
        """Test length of simple sequence."""
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        assert get_sequence_length(sequence) == 20

    def test_sequence_with_whitespace(self):
        """Test length of sequence with whitespace."""
        sequence = "ACDE FGH\nIKL MNP\tQRS TVW Y"
        assert get_sequence_length(sequence) == 20

    def test_empty_sequence(self):
        """Test length of empty sequence."""
        assert get_sequence_length("") == 0
        assert get_sequence_length("   \n\t  ") == 0


class TestCleanSequence:
    """Test class for clean_sequence function."""

    def test_clean_simple_sequence(self):
        """Test cleaning simple sequence."""
        sequence = "acdefghiklmnpqrstvwy"
        result = clean_sequence(sequence)
        assert result == "ACDEFGHIKLMNPQRSTVWY"

    def test_clean_sequence_with_whitespace(self):
        """Test cleaning sequence with whitespace."""
        sequence = "acde fgh\nikl mnp\tqrs tvw y"
        result = clean_sequence(sequence)
        assert result == "ACDEFGHIKLMNPQRSTVWY"

    def test_clean_empty_sequence(self):
        """Test cleaning empty sequence."""
        assert clean_sequence("") == ""
        assert clean_sequence("   \n\t  ") == ""

    def test_clean_mixed_case(self):
        """Test cleaning mixed case sequence."""
        sequence = "AcDeFgHiKlMn"
        result = clean_sequence(sequence)
        assert result == "ACDEFGHIKLMN"