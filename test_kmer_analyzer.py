import pytest
from kmer_analyzer import (
    validate_sequence,
    update_kmer_count,
    count_kmers_with_context,
    merge_kmer_data,
    write_results_to_file
)

# =========================
# Tests for validate_sequence
# =========================

def test_validate_sequence_valid():
    # Valid DNA sequence with proper characters and length
    assert validate_sequence("ATGTCTGAA", 2) is True


def test_validate_sequence_too_short():
    # Sequence shorter than k cannot form k-mers
    assert validate_sequence("AT", 3) is False


def test_validate_sequence_invalid_chars():
    # Sequence containing invalid characters should fail
    assert validate_sequence("ATG1TC", 2) is False


def test_validate_sequence_lowercase():
    # Lowercase letters are considered invalid
    assert validate_sequence("atgctg", 2) is False


def test_validate_sequence_empty():
    # Empty sequence is invalid
    assert validate_sequence("", 2) is False


def test_validate_sequence_exact_length():
    # Length == k gives no next character → invalid for this assignment
    assert validate_sequence("AT", 2) is False


# =========================
# Tests for update_kmer_count
# =========================

def test_update_kmer_count_new_kmer():
    # First time seeing a k-mer → count should be 1
    data = {}
    update_kmer_count(data, "AT", "G")
    assert data["AT"]["count"] == 1


def test_update_kmer_count_existing_kmer():
    # Repeated update should increment count
    data = {}
    update_kmer_count(data, "AT", "G")
    update_kmer_count(data, "AT", "G")
    assert data["AT"]["count"] == 2


def test_update_kmer_count_multiple_next_chars():
    # Same k-mer followed by different characters should track both
    data = {}
    update_kmer_count(data, "AT", "G")
    update_kmer_count(data, "AT", "C")

    assert data["AT"]["count"] == 2
    assert data["AT"]["next_chars"]["G"] == 1
    assert data["AT"]["next_chars"]["C"] == 1


# =========================
# Tests for count_kmers_with_context
# =========================

def test_count_kmers_with_context_basic():
    result = count_kmers_with_context("ATGTCTGAA", 2)

    # Ensure expected k-mers exist
    assert "AT" in result
    assert "TG" in result


def test_count_kmers_with_context_counts():
    result = count_kmers_with_context("ATGTCTGAA", 2)

    # TG appears twice
    assert result["TG"]["count"] == 2


def test_count_kmers_with_context_next_chars():
    result = count_kmers_with_context("ATGTCTGAA", 2)

    # AT is followed by G
    assert result["AT"]["next_chars"]["G"] == 1


def test_count_kmers_with_context_single_case():
    # Only one valid k-mer with context
    result = count_kmers_with_context("ATG", 2)

    assert result["AT"]["count"] == 1
    assert result["AT"]["next_chars"]["G"] == 1


# =========================
# Tests for merge_kmer_data
# =========================

def test_merge_kmer_data_basic():
    # Merging two datasets should combine counts
    a = {"AT": {"count": 1, "next_chars": {"G": 1}}}
    b = {"AT": {"count": 1, "next_chars": {"C": 1}}}

    merged = merge_kmer_data(a, b)

    assert merged["AT"]["count"] == 2
    assert merged["AT"]["next_chars"]["G"] == 1
    assert merged["AT"]["next_chars"]["C"] == 1


def test_merge_kmer_data_new_kmer():
    # New k-mer should be added
    a = {}
    b = {"TG": {"count": 2, "next_chars": {"A": 2}}}

    merged = merge_kmer_data(a, b)

    assert "TG" in merged
    assert merged["TG"]["count"] == 2


# =========================
# Tests for write_results_to_file
# =========================

def test_write_results_to_file(tmp_path):
    # Single k-mer output test
    kmer_data = {"AT": {"count": 1, "next_chars": {"G": 1}}}

    output = tmp_path / "output.txt"
    write_results_to_file(kmer_data, str(output))

    content = output.read_text()

    assert "AT 1 G:1" in content


def test_write_results_multiple_kmers(tmp_path):
    # Multiple k-mers should all appear formatted correctly
    kmer_data = {
        "AT": {"count": 2, "next_chars": {"G": 2}},
        "TG": {"count": 1, "next_chars": {"C": 1}}
    }

    output = tmp_path / "output.txt"
    write_results_to_file(kmer_data, str(output))

    content = output.read_text()

    assert "AT 2 G:2" in content
    assert "TG 1 C:1" in content


# =========================
# Full pipeline test 
# =========================

def test_full_pipeline(tmp_path):
    # Full real example from assignment
    sequence = "ATGTCTGTCTGAA"
    k = 2

    seq_data = count_kmers_with_context(sequence, k)
    global_data = merge_kmer_data({}, seq_data)

    output = tmp_path / "output.txt"
    write_results_to_file(global_data, str(output))

    result = sorted(output.read_text().splitlines())

    expected = sorted([
        "AT 1 G:1",
        "CT 2 G:2",
        "GA 1 A:1",
        "GT 2 C:2",
        "TC 2 T:2",
        "TG 3 A:1 T:2"
    ])

    assert result == expected
