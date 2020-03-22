import os
import json

from Bio import Entrez, SeqIO


def get_list_of_accession_numbers(logger, folder_path):
    file_paths = []
    with open(os.path.join(os.getcwd(), "sample_data.txt"), "r") as sample_file:
        sample_data = json.loads(sample_file.read())
        logger.log(f"Sequences found: {sample_data}")
        # Go through each record in the json file provided by the user
        for key, value in sample_data.items():
            Entrez.email = "k.patel1098@gmail.com"
            handle = Entrez.efetch(db="nucleotide", id=value, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            gene_file_path = os.path.join(folder_path, f"{key}.fasta")
            file_paths.append(gene_file_path)
            logger.log(f"Writing {key} fasta file")
            with open(gene_file_path, "w") as cdna_file:
                for feature in record.features:
                    # write each cds feature to a fasta file per gene
                    if feature.type == "CDS":
                        cdna_file.write(f">{feature.qualifiers['protein_id'][0]}\n")
                        cdna_file.write(f"{feature.location.extract(record).seq}\n")
    return file_paths


