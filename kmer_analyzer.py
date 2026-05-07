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
    
    # Check that sequence is long enough to form at least one k-mer
    if len(sequence) < k:
        return False
    
    # Ensure every charter is a valid DNA nucleotide
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
    
    # If k-mer has not been seen befor, initialize its data structure
    if kmer not in kmer_data:
        kmer_data[kmer] = {'count': 0, 'next_chars': {}}
    
    # Increment total count for this k-mer
    kmer_data[kmer]['count'] += 1
    
    # If this next charater hasn't been seen for this k-mer, initialize it
    if next_char not in kmer_data[kmer]['next_chars']:
        kmer_data[kmer]['next_chars'][next_char] = 0
    
    # Increment the frequency of this next character
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
    
    # Dictionary to store reults for this single sequence
    kmer_data = {}
    
    # Loop through sequence to extract all k-mers with a following character
    # Stop at len (sequence) - k so it does not go out of bounds
    for i in range(len(sequence) - k):
        
        # Extract k-mer substring
        kmer = sequence[i:i+k]
        
        # Get the character immediately following the k-mer
        next_char = sequence[i+k]
        
        # Update counts using helper function
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
    
    # Loop through each k-mer in the new dataset
    for kmer in new_data:
        
        # If k-mer doesn't exist in global data, initialize it
        if kmer not in global_data:
            global_data[kmer] = {'count': 0, 'next_chars': {}}

        # Add total k-mer count from this sequence 
        global_data[kmer]['count'] += new_data[kmer]['count']

        # Merge next-character frequencies 
        for char, freq in new_data[kmer]['next_chars'].items():
            
            # Initialize character if not already present
            if char not in global_data[kmer]['next_chars']:
                global_data[kmer]['next_chars'][char] = 0

            # Add frequency from this sequence 
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
   
    # Sort k-mers alphabetically for consistent output
    sorted_kmers = sorted(kmer_data.keys())
    
    # Open output file for writing
    with open(output_filename, 'w') as f:
        
        # Loop through each k-mer
        for kmer in sorted_kmers:
            
            # Total count of this k-mer
            total = kmer_data[kmer]['count']
            
            # Dictionary of next-character frequencies
            next_chars = kmer_data[kmer]['next_chars']
            
            # Convert next-character dictionary into formatted string
            # Example: "A:1 T:2"
            next_char_str = " ".join(
                f"{char}:{freq}" 
                for char, freq in sorted(next_chars.items())
            )
            
            # Write formmated line to file
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
   
    # Read command line arguments
    sequence_file = sys.argv[1]
    k = int(sys.argv[2])
    output_file = sys.argv[3]
    
    print(f"Reading sequences from {sequence_file}...")
    
    # Dictionary to hold combined results across all sequences
    global_kmer_data = {}

    # Open and read the input file line by line
    with open(sequence_file, 'r') as f:
        for sequence in f:
            
            # Remove whitespace/newline characters
            sequence = sequence.strip()

            # Skip invalid sequences and print warning
            if not validate_sequence(sequence, k):
                print(f"  Warning: Skipping sequence")
                continue
            
            # Process individual sequence to get k-mer counts
            seq_data = count_kmers_with_context(sequence, k)
            
            # Merge this sequence's data into the global dataset
            global_kmer_data = merge_kmer_data(global_kmer_data, seq_data)
            
        # After all sequences are processed, write final results once
        write_results_to_file(global_kmer_data, output_file)


# Entry point: ensures main() runs only when script is executed directly 
if __name__ == '__main__':
    main()
