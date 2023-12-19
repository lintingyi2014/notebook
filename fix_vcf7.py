#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import sys
import pysam

class Application:
    def __init__(self):
        self.args = Application.get_args()
        self.vcf = pysam.VariantFile(sys.stdin)
        self.fasta = pysam.FastaFile(self.args.fasta)

    @staticmethod
    def get_args() -> argparse.Namespace:
        """
        Read command line arguments and set input and output sources.

        :return: parsed command line arguments
        """

        parser = argparse.ArgumentParser(prog="fixvcf.py", description="Fixes deletion represented by -")
        parser.add_argument('--version', action='version', version='%(prog)s 0.1')

        parser.add_argument("input", help="input file, use - to read from stdin")
        parser.add_argument("-f", "--fasta", help="fasta reference file", required=True)
        parser.add_argument("-o", "--output", help="output file")

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

        args = parser.parse_args()

        if args.output:
            sys.stdout = open(args.output, "w")

        if args.input != "-":
            sys.stdin = open(args.input, "r")

        return args

    def fix_record(self, record: pysam.libcbcf.VariantRecord) -> pysam.libcbcf.VariantRecord:
        if record.alts is not None and "-" in record.alts:
            ref = self.fasta.fetch(region=f"{record.chrom}:{record.pos - 1}-{record.pos - 1}")    
            record.pos = record.pos - 1
            record.ref = f"{ref}{record.ref}"
            record.alts = tuple(map(lambda alt: ref if alt == "-" else f"{ref}{alt}", record.alts))
        
        if record.ref is not None and "-" in record.ref:
            ref = self.fasta.fetch(reference=record.chrom, start=record.pos - 1, end=record.pos)

            record.pos = record.pos
            record.ref = f"{ref}"
            record.alts = tuple(map(lambda alt: f"{ref}{alt}", record.alts))
        return record

    def start(self):
        print(self.vcf.header, end="")
        for record in self.vcf:
            if record.alts is not None:  # Check if alts is not None
                fixed_record = self.fix_record(record)
                print(fixed_record, end="")

if __name__ == "__main__":
    app = Application()
    app.start()

