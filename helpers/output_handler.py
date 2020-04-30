import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation


def output_transcript_as_fasta(global_data, transcript):
    global_data.logger.log(f"Writing transcript cDNA fasta file for {transcript.transcript_name}")
    transcript_dir_path = os.path.join(global_data.folder_path, f"{transcript.transcript_name}")
    try:
        os.mkdir(transcript_dir_path)
    except Exception as e:
        pass
    cdna_path = os.path.join(transcript_dir_path, f"{transcript.transcript_name}_transcript.fasta")
    with open(cdna_path, "w") as cdna_file:
        cdna_file.write(f">{transcript.transcript_id}\n")
        cdna_file.write(f"{transcript.coding_sequence}\n")


def _handle_codons_in_features(input_sequence):
    if input_sequence is not None:
        remainder = len(input_sequence) % 3

        if remainder == 1:
            input_sequence += "GG"
        elif remainder == 2:
            input_sequence += "G"
        else:
            pass
    return input_sequence


def export_genbank_files(global_data):
    for transcript in global_data.ensembl_gene.transcripts:
        if transcript.transcript_name in global_data.valid_transcripts:
            global_data.logger.log(f"Exporting {transcript.transcript_name} as genbank")
            base_count = 0

            promoter_seq = _handle_codons_in_features(str(global_data.promoter.seq))
            five_utr_seq = _handle_codons_in_features(str(transcript.five_prime_utr_sequence))
            final_seq = ""

            promoter_feature = None
            if global_data.promoter is not None and len(promoter_seq) > 0:
                promoter_length = len(promoter_seq)
                final_seq += str(promoter_seq)
                promoter_start = base_count
                promoter_feature = SeqFeature(FeatureLocation(start=promoter_start, end=base_count + promoter_length),
                                              type=f"HSP70 promoter - {global_data.promoter.id}")
                base_count += promoter_length

            five_utr_feature = None
            if five_utr_seq is not None and len(five_utr_seq) > 0:
                utr_length = len(five_utr_seq)
                final_seq += str(five_utr_seq)
                utr_start = base_count
                five_utr_feature = SeqFeature(FeatureLocation(start=utr_start, end=base_count + utr_length),
                                              type="5' UTR sequence")
                base_count += utr_length

            coding_seq_start = base_count
            base_count += len(transcript.coding_sequence)
            final_seq += str(transcript.coding_sequence)

            coding_sequence_feature = SeqFeature(FeatureLocation(start=coding_seq_start, end=base_count),
                                                 type=f"{transcript.transcript_name} Coding Sequence")

            coding_sequence = Seq(final_seq, IUPAC.unambiguous_dna)
            record = SeqRecord(coding_sequence,
                               id='123456789',
                               name=f'{transcript.transcript_name} Gene Fusion'.replace(" ", "_"),
                               description=f'{transcript.transcript_name} fusion with promoter'.replace(" ", "_"))
            if promoter_feature is not None:
                record.features.append(promoter_feature)
            if five_utr_feature is not None:
                record.features.append(five_utr_feature)
            if coding_sequence_feature is not None:
                record.features.append(coding_sequence_feature)

            output_file = open(os.path.join(global_data.folder_path, f"{transcript.transcript_name}",
                                            f'{transcript.transcript_name}.gb'), 'w')
            output_fasta_file = open(os.path.join(global_data.folder_path, f"{transcript.transcript_name}",
                                                  f'{transcript.transcript_name}_fusion.fasta'), 'w')
            SeqIO.write(record, output_file, 'genbank')
            SeqIO.write(SeqRecord(coding_sequence, f"{transcript.transcript_name} - {transcript.transcript_id}",
                                  description="Fusion"), output_fasta_file, 'fasta')
