## K-mer Analyzer using Python
Analyze k-mers in genome sequence fragments using Python

## Running the Script
python kmer_analyzer.py <sequence_file><output_file>

## Output Format
The output file contains one line per k-mer, sorted alphabetically. Each line has the following format:
<kmer> <total_count> <char1>:<count1> <char2>:<count2>

## Data Dictionary
- kmer: The substring of length k
- total_count: The number of times the k-,er appears in all sequences
- char:count pairs: The frequency of each nucleotide that immediately follows the k-mer

## AI Usage Statement 
This project used AI assistance (Copilot) to help debug code, imprvode structure, and generate test cases. All code was reviewed, tested, and understood be author.

