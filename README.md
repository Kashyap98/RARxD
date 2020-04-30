# Gene Fusion Tools - RARxD

## Purpose - Project Abstract
Gene fusion is a useful method for researchers to study the expression and functions of various genes in eukaryotic organisms. Researchers that study zebrafish (Danio rerio) use the HSP70 promoter to study the expression of genes through temperature control. Promoter-gene fusion allows scientists to fuse genes that are not naturally influenced by the HSP70 promoter to the target gene. In order to create a proper fusion, one must determine which gene transcript they would like to fuse to the promoter region. An open source software tool can be useful to aggregate data necessary to make this decision. This gene fusion tool, available at https://github.com/Kashyap98/gene_fusion_tool, aggregates Ensembl transcript annotations and gene identity in NCBI’s BLAST databases to assist the scientist in creating the proper fusion. The tool outputs this information in a standardized format for comparison as well as a properly fused sequence so the researcher may perform Gibson Assembly with the necessary primers.

## Usage
#### There are 2 main ways to use the Gene Fusion Tool (Command Line, or GUI)
#### Command Line
| Argument | Description |
| -------- | ----------- |
| --help | Show help message text |
| --name (required) |  a name for your test run |
| --gene_id (required)| Ensembl Accession ENSDARGXXXXXXXXXXX: Enter the id found in the ensembl database|
| --promoter_id (required)| Promoter accession from NCBI to prepend to the transcript. |
| --quiet  | y to only log to log file, n (default) to log to log file and print to console. |
| --blast_type | LOCAL (default) to local blastn database, REMOTE to remotely BLAST NCBI server. |

Example command: \
`python controller.py --name raraa --gene_id ENSDARG00000056783 --promoter_id AF158020.1` \
This command creates a new folder named run_raraa in the project directory. It is accessing the raraa gene (ENSDARG00000056783) from Ensembl and
using the HSP70 promoter from NCBI (AF158020.1). It will default to a LOCAL BLASTn and log output to a file and the console.

#### GUI

We were unable to create an installer in the scope of this project. Therefore to run the Gene Fusion Tool with the GUI the following command must be run. \
`python main_dialog.py` \
Entering raraa into the Run Name input, ENSDARG00000056783 in the Ensembl Accession input, and clicking start run will yield the same results as above.
Other features that can be changed are the REMOTE/LOCAL BLAST and changing the filters used. \
Once the run is complete, if it does not show up on the right click the refresh button. The view button will bring up the logs, input data/settings, and the valid transcripts and their results. There are also clickable links that can take you to your results.

## Installation

#### Python
This was made and tested in python 3.6. \
Python can be installed from here: \
https://www.python.org/downloads/release/python-368/

#### Retrieve Project
Using the command line, navigate to the folder you wish to add the project to. \
Type `git clone https://github.com/Kashyap98/gene_fusion_tool` to get project files. \
Type `cd gene_fusion_tool` to enter project directory.

From here, controller.py is the most important file for the operation of this project. For these examples,
if you have multiple installations of python use `python3 or pip3` instead of `python or pip`.

To install most of the dependencies run. `pip install -r requirements.txt`

#### NCBI BLAST+
NCBI's BLAST+ must be installed separately (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) \
Select and install the appropriate version for your operating system.

#### BioPython
Installed using pip (requirements.txt) (https://biopython.org/) \
BioPython Manual \
http://biopython.org/DIST/docs/tutorial/Tutorial.html

#### Local BLASTn Database
Download this version of the zebrafish genome
(https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Danio_rerio/all_assembly_versions/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_cds_from_genomic.fna.gz) \
Unzip the gz file and rename it to zebrafish_cdna_from_genomic.fasta \
Make sure you have BLAST+ installed from above. \
Run this command in the correct directory. \
`makeblastdb -in zebrafish_cdna_from_genomic.fasta -parse_seqids -title zebrafish_blastdb -dbtype nucl`

#### Future Directions:
The future steps of developing this application are numerous. One can modify the inputs to allow the program to be more generalized and not just specialize in Danio rerio. In addition, one can modify the inputs to allow for multiple runs to be done at once rather than one at time. Furthermore, more filters can be added to further specify the different qualities that make a transcript undesirable. Other paths for development include adding other types of BLAST, such as BLASTx and allowing the user to decide if one, or both, should be run. In addition, scripts can be written to make the installation of this GUI program, pyensembl package, and BLAST database easier. Furthermore, work can be done in optimizing what, and how, the generated data is displayed to the user. Additionally, another goal may be to optimize the communication the Gene Fusion tool has with the user while it is working in GUI mode rather than the command line. 