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
    -d [input] : Input file.
    -w [working folder] : Directory for temporary files.
    -t [# of threads] : Number of CPU threads.
    -o [output] : Output file name.
    -n [# of candidates] : Number of candidates to consider for gapped extension.
    -a [overlap size] : Minimum overlap size.
    -g [0/1] : Whether to output gapped extension start points.
    -x [0/1] : Sequencing platform (0 for PacBio, 1 for Nanopore).

Note: The -n parameter is crucial for optimizing FODI's performance. It's recommended to adjust this parameter based on genome size and read length.
## Demo

This section provides a quick demonstration of how to use our long-read overlap detection tool with a sample dataset.

## Sample Data

For this demonstration, we will use a subset of PacBio reads from an *E. coli* dataset. You can download the FASTQ file from https://goo.gl/Z75V5R.

## Running the Tool

To run the tool, follow these steps:

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/EghdamiKerg/FODI.git
    ```

2.  **Navigate into the project directory.**
    (Assuming `FODI` is the project directory name after cloning)
    ```bash
    cd FODI
    ```

3.  **Open the Code::Blocks project file:**
    Launch Code::Blocks IDE and open the `.cbp` file located in the root of the project directory.

4.  **Compile the project:**
    In Code::Blocks, go to `Build` > `Build` or press `Ctrl+F9`. This will compile the source code and generate the executable.

5.  **Navigate to the executable directory:**
    Open your terminal and change directory to the location of the compiled executable:
    ```bash
    cd bin/Debug/
    ```

6.  **Run the executable with desired parameters:**
    Execute the program from the terminal, providing the necessary parameters.

## Example Usage

To run the overlap detection tool with the provided sample data, use the following parameters:

```bash
./firstCodeBlocksProject -j 0 -d /path/to/ecoli.fastq -o output.can -w owrk -x 0 -n 100 -a 2000 -t 8   ```

Where:

    -j 0: specifies the task is detecting overlapping candidates

    -d /path/to/ecoli_clipped.fastq: Specifies the path to the input FASTQ file. Replace /path/to/ecoli_clipped.fastq with the actual path to the downloaded ecoli_clipped.fastq file.

    -o output.can: Specifies the output file name (in can format).

    -w owrk: specifies the path to the Directory of temporary files

    -x 0: Pacbio sequencing platform is used

    -n 100: specifies number of candidates

    -a 2000: specifies Minimum overlap size

## Expected Output

If the job is detecting overlapping candidates, the results are output in can format, with each result on one line and containing 9 fields:

[A ID] [B ID] [A strand] [B strand] [A gapped start] [B gapped start] [voting score] [A length] [B length]

