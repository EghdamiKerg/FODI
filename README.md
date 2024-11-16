# FODI
FODI: Fast Overlap Detection with Dual Indexing

FODI is a fast and accurate tool for detecting overlaps between long reads. It utilizes a novel dual-indexing approach to efficiently identify candidate overlaps and a dynamic programming algorithm to refine the overlap regions.

Input:

    FASTQ format for input reads

Output:

    Overlap candidates in CAN format
    Pairwise alignments in M4 format

Options:

    -j [task] (0 or 1): Select between detecting overlapping candidates only (0) or outputting pairwise alignments in M4 format (1).
    -d [file] : Input file.
    -w [working folder] : Directory for temporary files.
    -t [# of threads] : Number of CPU threads.
    -o [output] : Output file name.
    -n [# of candidates] : Number of candidates to consider for gapped extension.
    -a [overlap size] : Minimum overlap size.
    -g [0/1] : Whether to output gapped extension start points.
    -x [0/1] : Sequencing platform (0 for PacBio, 1 for Nanopore).

Note: The -n parameter is crucial for optimizing FODI's performance. It's recommended to adjust this parameter based on genome size and read length.
