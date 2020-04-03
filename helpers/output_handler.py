import os


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
