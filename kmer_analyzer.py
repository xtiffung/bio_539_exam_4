import sys

def validate_sequence(sequence, k):
    """Validate a DNA sequence.

    Parameters:
        sequence (str): The DNA sequence string.
        k (int): Length of k-mers.

    Returns:
        bool: True if sequence is valid, False otherwise.

    A valid sequence:
    - Has length >= k
    - Contains only valid nucleotides: A, C, G, T
    """
    if len(sequence) < k:
        return False
    for nucleotide in sequence:
        if nucleotide not in 'ACGT':
            return False
    return True

def update_kmer_count(kmer_data, kmer, next_char):
    """Update k-mer counts and next-character frequencies.

    Parameters:
        kmer_data (dict): Dictionary storing k-mer information.
        kmer (str): Current k-mer.
        next_char (str): Character following the k-mer.

    Returns:
        dict: Updated kmer_data dictionary.

    This function:
    - Increments total k-mer count
    - Tracks how often each nucleotide follows the k-mer
    """
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}
    
    kmer_data[kmer]['count'] += 1
    
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0
    kmer_data[kmer]['next_chars'][next_char] += 1

    return kmer_data

def count_kmers_with_context(sequence, k):
    """
    Extract k-mers and record following nucleotide frequencies.

    Parameters:
        sequence (str): DNA sequence.
        k (int): Length of k-mers.

    Returns:
        dict: Dictionary containing k-mer counts and next-character data.

    Iterates through the sequence and:
    - Extracts each k-mer
    - Records the character immediately following it
    """
    kmer_data = {}
    
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        next_char = sequence[i+k]
        
        kmer_data = update_kmer_count(kmer_data, kmer, next_char)
    
    return kmer_data

def merge_kmer_data(global_data, new_data):
    """Merge k-mer statistics from one dataset into another.

    Parameters:
        global_data (dict): Aggregated k-mer data.
        new_data (dict): k-mer data from one sequence.

    Returns:
        dict: Updated global k-mer data.

    Combines:
    - Total k-mer counts
    - Next-character frequencies
    """
    for kmer in new_data:
        if kmer not in global_data:
            global_data[kmer] = {'count': 0, 'next_chars': {}}

        # Add total k-mer count
        global_data[kmer]['count'] += new_data[kmer]['count']

        # Add next character counts
        for char, freq in new_data[kmer]['next_chars'].items():
            if char not in global_data[kmer]['next_chars']:
                global_data[kmer]['next_chars'][char] = 0

            global_data[kmer]['next_chars'][char] += freq

    return global_data

def write_results_to_file(kmer_data, output_filename):
    """Write k-mer statistics to an output file.

    Parameters:
        kmer_data (dict): Dictionary of k-mer data.
        output_filename (str): Output file path.

    Output format:
        <kmer> <total_count> <char1>:<count1> <char2>:<count2> ...

    K-mers are sorted alphabetically for consistent output.
    """
    sorted_kmers = sorted(kmer_data.keys())
    
    with open(output_filename, 'w') as f:
        for kmer in sorted_kmers:
            total = kmer_data[kmer]['count']
            next_chars = kmer_data[kmer]['next_chars']
            
            next_char_str = " ".join(
                f"{char}:{freq}" 
                for char, freq in sorted(next_chars.items())
            )
            
            f.write(f"{kmer} {total} {next_char_str}\n")

def main():
    """Main function to execute the k-mer analysis pipeline.

    Steps:
    1. Parse command-line arguments
    2. Read input sequences from file
    3. Validate sequences
    4. Compute k-mer statistics
    5. Merge results across sequences
    6. Write final output to file
    """
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]
    
    print(f"Reading sequences from {sequence_file}...")
    
    # create global storage
    global_kmer_data = {}

    with open(sequence_file, 'r') as f:
        for sequence in f:
            sequence = sequence.strip()

            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence")
                continue
            
            # process one sequence
            seq_data = count_kmers_with_context(sequence, k)
            
            # merge results 
            global_kmer_data = merge_kmer_data(global_kmer_data, seq_data)
            
        # write once at the end
        write_results_to_file(global_kmer_data, output_file)


if __name__ == '__main__':
    main()
