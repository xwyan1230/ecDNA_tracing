"""
Performs some quality filtering on FASTQ reads.
"""

import argparse
# from lib2to3.pytree import convert
import numpy as np
import pysam


def mask_low_quality(fastq_fn: str, output_fn: str, min_quality: int = 25, verbose: bool = False):

    _iter = 0

    with pysam.FastxFile(fastq_fn) as parental_fastq, open(
        output_fn, mode="w"
    ) as masked_fastq:
        for record in parental_fastq:

            sequence, quality = (
                np.array([char for char in record.sequence]),
                np.array(
                    [convert_character_to_quality(char) for char in record.quality]
                ),
            )
            sequence[(quality < min_quality)] = "N"
            record.sequence = "".join(sequence)
            masked_fastq.write(str(record) + "\n")

            _iter += 1
            if verbose and _iter % 1e6 == 0:
                print(f"Processed {_iter} records.")


def convert_character_to_quality(character: str) -> int:
    """Converts a ASCII character from FASTQ to Q score."""
    ascii = ord(character)
    return ascii - 33


def main():

    parser = argparse.ArgumentParser(description="Quality filter reads.")
    parser.add_argument("r1", type=str, help="Path to R1 FASTQ")
    parser.add_argument("r2", type=str, help="Path to R2 FASTQ")
    parser.add_argument("--min_quality", type=int, default=25)
    parser.add_argument("--verbose", action='store_true')

    args = parser.parse_args()

    r1_path = args.r1
    r2_path = args.r2
    min_quality = args.min_quality
    verbose = args.verbose

    r1_masked = r1_path.replace(".fastq.gz", ".masked.fastq")
    r2_masked = r2_path.replace(".fastq.gz", ".masked.fastq")

    mask_low_quality(r1_path, r1_masked, min_quality=min_quality, verbose=verbose)
    mask_low_quality(r2_path, r2_masked, min_quality=min_quality, verbose=verbose)


if __name__ == "__main__":
    main()