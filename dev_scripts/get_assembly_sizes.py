#!/usr/bin/env python3
import os
import json
import csv

def extract_assembly_info(jsonl_file):
    with open(jsonl_file, 'r') as file:
        data = json.load(file)

        accession = data.get("accession")
        total_sequence_length = data.get("assemblyStats", {}).get("totalSequenceLength")

        return accession, total_sequence_length

def process_directory(directory):
    results = []

    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".jsonl"):
                file_path = os.path.join(root, file)
                accession, total_sequence_length = extract_assembly_info(file_path)
                results.append((accession, total_sequence_length))

    return results

def write_to_tsv(output_file, data):
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(["Assembly Accession", "Total Sequence Length"])
        writer.writerows(data)

# Example usage
directory_to_search = "/scratch/molevo/dts/ODP_genomes/GenDB_unannotated/odp_ncbi_genome_db/output/source_data/unannotated_genomes/"
output_tsv_file = "output_assembly_sizes.tsv"

results = process_directory(directory_to_search)
write_to_tsv(output_tsv_file, results)
