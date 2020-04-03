import jsonpickle
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML, NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
import os


# organize the blast result in an easily accessible class
class BlastResult:

    def __init__(self, global_data, result_string, transcript, record=None, is_local=True, is_writing=True):
        if is_local:
            accession, match_length, bitscore = self.parse_result_string(result_string)
        else:
            accession, match_length, bitscore = result_string.split(",")
        self.accession = accession
        self.match_length = match_length
        self.bitscore = bitscore
        self.transcript = transcript
        self.record = record
        self.global_data = global_data
        if is_writing:
            self.export_as_json()

    def parse_result_string(self, input_string):
        single_result_list = input_string.split(',')
        match_length = single_result_list[1]
        bitscore = single_result_list[2]
        accession = single_result_list[0].split("_cds_")[1].split(".")[0]
        return accession, match_length, bitscore

    def export_as_json(self):
        with open(os.path.join(self.global_data.folder_path, f"{self.transcript.transcript_name}",
                               f"{self.accession}.json"), "w") as output_file:
            output_file.write(jsonpickle.encode(self))


# Perform the blast
def blast(global_data, transcript):
    transcript_path = os.path.join(global_data.folder_path, f"{transcript.transcript_name}",
                                   f"{transcript.transcript_name}_transcript.fasta")
    blastn = NcbiblastnCommandline(query=transcript_path, db=os.path.join(os.getcwd(),
                                                                          "zfish_cdna_from_genomic.fna"),
                                   outfmt='"10 sseqid length bitscore"',
                                   max_target_seqs=5)
    result = list(blastn())
    result_list = result[0].split("\n")
    blast_results = []

    if len(result_list) > 0:
        for result in result_list:
            if result != "":
                blast_result = BlastResult(global_data, result, transcript, is_local=True, is_writing=True)
                blast_results.append(blast_result)
    else:
        blast_results = None

    return blast_results


def remote_blast(global_data, transcript):
    global_data.logger.log(f"BLASTing {transcript.transcript_name} online... This may take a while")
    blast_results = []
    result_handle = NCBIWWW.qblast("blastn", "nr", transcript.coding_sequence, megablast=True,
                                   entrez_query="Danio[genus]")
    blast_records = list(NCBIXML.parse(result_handle))
    if len(blast_records) > 0:
        result = blast_records[0]
        for alignment in result.alignments:
            result_string = f"{alignment.accession},{alignment.length},{alignment.hsps[0].bits}"
            blast_result = BlastResult(global_data, result_string, transcript, alignment, is_local=False,
                                       is_writing=True)
            blast_results.append(blast_result)
    else:
        blast_results = None

    return blast_results


def handle_transcripts(global_data, is_remote=False):
    transcript_results = []
    for transcript in global_data.ensembl_gene.transcripts:
        if transcript.transcript_name not in global_data.valid_transcripts:
            continue
        if is_remote:
            transcript_result = remote_blast(global_data, transcript)
        else:
            transcript_result = blast(global_data, transcript)
        if transcript_result is not None:
            transcript_results.append(transcript_result)
    return transcript_results


def locally_handle_blasting_transcripts(global_data):
    global_data.logger.log("Working on BLASTing transcripts locally")
    transcript_results = handle_transcripts(global_data, is_remote=False)
    global_data.transcript_results = transcript_results
    return global_data


def remotely_handle_blasting_transcripts(global_data):
    global_data.logger.log("Working on BLASTing transcripts remotely")
    transcript_results = handle_transcripts(global_data, is_remote=True)
    global_data.transcript_results = transcript_results
    return global_data


def get_blasted_transcript_information(global_data):
    for transcript_list in global_data.transcript_results:
        for transcript in transcript_list:
            Entrez.email = "k.patel1098@gmail.com"
            handle = Entrez.efetch(db="protein", id=transcript.accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            transcript.record = record
    return global_data
