#!/usr/bin/env python3
"""
Author: George David Apostolidis
Student ID: 1746227
Description: Create a multiple sequence alignment of the protein sequences using
            the program Clustal Omega
Usage: Usage:python3 exam.py ins_mammalia.fasta INS_BOVIN
"""
import subprocess
import sys
from typing import List, Dict
from dataclasses import dataclass



@dataclass
class AlignedSequence:
    swissprot_seq: str
    entry_id: str
    entry_name: str
    bases: str

@dataclass
class ClustalOutput:
    sequences: List[AlignedSequence]


def input_args():
    """
        Parse and validate command-line arguments.

        Returns
        -------
        filenames : List[str]
            List containing the input FASTA filename.
        reference_entry: str
            The reference entry given by the user
    """

    if len(sys.argv) != 3:
        print("Usage:python3 exam.py ins_mammalia.fasta INS_BOVIN")
        sys.exit(1)

    filename = sys.argv[1]
    reference_entry = sys.argv[2]
    return filename, reference_entry


def run_cmd(cmd: List[str], filename)-> None:
    """
        Run a command using subprocess and handle errors.

        Parameters
        ----------
        filename: Name of file to run.
        cmd : List[str]
            Command and its arguments to be executed.

        Returns
        -------
        subprocess.CompletedProcess
            The result of the executed command.
    """

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd)}", file=sys.stderr)
        print(f"Return code: {e.returncode}", file=sys.stderr)
        print(f"Output: {e.output}", file=sys.stderr)
        sys.exit(1)

def parse_clustal_output(OUTPUT_FILE, clustal_output=None):
    """
        Parse the Clustal output file to extract alignment sequences.

        Parameters
        ----------
        OUTPUT_FILE : str
            Path to the Clustal output file.
        clustal_output : ClustalOutput, optional
            Existing ClustalOutput object to populate (default is None).

        Returns
        -------
        ClustalOutput
            Parsed ClustalOutput object containing alignment sequences.
    """
    with open(OUTPUT_FILE, "r") as f:
        lines = f.readlines()
        count = 0
        sequences = []
        for i, line in enumerate(lines):
            line = line.strip()

            if line.startswith(">"):
                aligned_seq = line.split("|")

                swissProt_db_name = aligned_seq[0].replace(">","")
                entry_id = aligned_seq[1]
                entry_name = aligned_seq[2].split(" ")[0]
                bases = ""

            elif "CLUSTAL" not in line and "*" not in line and "." not in line and i < len(lines)-1:
                bases += line.rstrip("\n")
                next_Line = lines[i+1]
                if len(next_Line) > 1 and not next_Line.startswith(">"):
                    bases+=next_Line.rstrip("\n")

            count += 1
            sequence = AlignedSequence(
                swissprot_seq=swissProt_db_name,
                entry_id=entry_id,
                entry_name=entry_name,
                bases=bases
            )
            if count % 3 == 0:
                sequences.append(sequence)

    clustal_output = ClustalOutput(sequences=sequences)

    return clustal_output


def find_conserved_positions(clustal_output):
    """
        Finds conserved positions in the alignment where all sequences have the same amino acid

        Parameters
        ----------
        clustal_output : ClustalOutput, optional
            Existing ClustalOutput object to populate (default is None).

        Returns
        -------
        same_amino : A dictionary with the indexes of the same amino acids
    """
    same_amino = {}
    for sequence in clustal_output.sequences:
        for i,base in enumerate(sequence.bases):
            for j, base2 in enumerate(sequence.bases):
                if any(char in base for char in base2) and i!=j and not base2 in same_amino.keys():
                    same_amino[base2] = i

    with open("out.txt", "w", encoding="utf-8") as fh:
        for key in same_amino:
            string_to_write = str(same_amino[key]) + "-->" + key + "\n"
            fh.write(str(string_to_write))

    return same_amino

def find_aligned_sequence(reference_entry, clustal_output,same_amino_acids):
    """
        Finds aligned sequence with a reference entry

        Parameters
        ----------
        same_amino_acids :  A dictionary with the indexes of the same amino acids
        reference_entry : The reference entry to check the aligned sequence
        clustal_output : ClustalOutput, optional
            Existing ClustalOutput object to populate (default is None).
    """
    found = False
    for sequence in clustal_output.sequences:
        if reference_entry == sequence.entry_name:
            found = True
            base_to_be_printed = sequence.bases.replace("-","").lower()

            for i,char in enumerate(base_to_be_printed):
                for key in same_amino_acids:
                    if char==key.lower():
                        base_to_be_printed = base_to_be_printed.replace(char,key)

            print(reference_entry,"-->",base_to_be_printed)

    if not found:
        print(reference_entry,"--> Entry not found",)


def main():
    """
        Main function to execute the Clustalo Command and process results.
    """

    filename, reference_entry = input_args()

    cmd = [
        "clustalo",
        "-i", filename,
        "-o",filename.split(".")[0] + ".clustal",
        "--force"
    ]

    run_cmd(cmd, filename)

    clustal_output = parse_clustal_output(filename.split(".")[0] + ".clustal")
    same_amino_acids = find_conserved_positions(clustal_output)
    find_aligned_sequence(reference_entry,clustal_output,same_amino_acids)


if __name__ == "__main__":
    main()

