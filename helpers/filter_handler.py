from helpers.output_handler import output_transcript_as_fasta


def __remove_transcript_helper(global_data, transcript, removal_message):
    global_data.logger.log(f"Removing {transcript.transcript_name} for reason: {removal_message}")
    global_data.invalidate_transcript(transcript.transcript_name)


def _validity_checker_boolean(transcript_variable):
    is_valid = True
    if transcript_variable is False:
        is_valid = False
    return is_valid


def filter_transcripts(global_data):
    global_data.logger.log("Working on filtering transcripts")
    for transcript in global_data.ensembl_gene.transcripts:
        global_data.logger.log(f"Working on filtering transcript {transcript.transcript_name}")

        # validate chosen parameters, later this will be done using settings
        if _validity_checker_boolean(transcript.complete) is False:
            __remove_transcript_helper(global_data, transcript, "Transcript is not complete")
            continue

        if _validity_checker_boolean(transcript.is_protein_coding) is False:
            __remove_transcript_helper(global_data, transcript, "Transcript is not protein coding")
            continue

        if _validity_checker_boolean(transcript.contains_start_codon) is False:
            __remove_transcript_helper(global_data, transcript, "Transcript does not contain start codon")
            continue

        if _validity_checker_boolean(transcript.contains_stop_codon) is False:
            __remove_transcript_helper(global_data, transcript, "Transcript does not contain stop codon")
            continue

        output_transcript_as_fasta(global_data, transcript)

    return global_data
