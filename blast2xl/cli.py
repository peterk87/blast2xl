# coding: utf-8
"""Console script for blast2xl."""
import logging
import re
import sys
from pathlib import Path
from typing import Dict, Tuple, Set

import click

from blast2xl.blast2xl import blast_tsv_to_df, output_xlsx_report
from blast2xl.io import write_seqs_to_taxonomy_dirs
from blast2xl.log import init_console_logger

logger = logging.getLogger(__name__)


@click.command()
@click.option('--blast-tsv-dir', type=click.Path(exists=True), required=True,
              help='Input directory with BLAST tab-delimited (TSV) files. The base filename should be the sample '
                   'name. If a `--seq-dir` is provided, the BLAST TSV base filenames should match the FASTA base '
                   'filenames.')
@click.option('--blast-tsv-sample-name-pattern', type=str, default=None,
              help='Regex pattern to extract sample name from BLAST TSV filename. '
                   'For example "^blastn-(.+)-vs-nt\\.tsv" to match "blastn-whatever-you-want-123-vs-nt.tsv" to '
                   'pull out sample name "whatever-you-want-123"')
@click.option('--seq-dir', type=click.Path(exists=True),
              help='Input directory with FASTA sequences from BLAST results. The base filename should be the sample '
                   'name and match the base filename for each BLAST result file.')
@click.option('--top-n-results', default=-1, type=int,
              help='Only take the top N results for each sequence when generating the "All BLAST Results" sheet in '
                   'the XLSX report.')
@click.option('-o', '--output-xlsx', type=click.Path(),
              required=True,
              help='XLSX output file path')
@click.option('-O', '--seq-outdir', type=click.Path(), default='taxonomy-sorted-sequences',
              help='Taxonomy organized sequence output directory')
@click.option('--keep-orientation', is_flag=True,
              help='Keep original orientation of BLAST searched sequences when outputting to `--seq-outdir`. Default '
                   'is to save in "plus" strand orientation to keep sequence orientations consistent with most NCBI '
                   'deposited sequences.')
@click.option('-v', '--verbose', default=0, count=True, help='Logging verbosity')
def main(blast_tsv_dir,
         blast_tsv_sample_name_pattern,
         seq_dir,
         top_n_results,
         output_xlsx,
         seq_outdir,
         keep_orientation,
         verbose):
    """blast2xl: BLAST XLSX Report Creator

    Given tabular (outfmt=6) blastn output with a few extra fields on top of the
    default standard ('std') set:

    \b
    \033[2m$\033[0m blastn ... -outfmt "6 qaccver saccver pident length mismatch gapopen qstart \
    qend sstart send evalue bitscore \033[32mqlen slen stitle staxid ssciname\033[0m" ...

    This script will generate an Excel XLSX report with the following sheets:

    \b
    - Sheet 1 "Top BLAST Results":
        - top hit by bitscore for each sequence in each input file only sorted by
          filename basename without extension and by sequence ID
    - Sheet 2 "All BLAST Results":
        - all hits sorted by filename, sequence ID and bitscore

    Optionally, this script can also output plaintext tab-delimited files for each
    sheet if specified by the user.

    Sequences from the different samples can also be output into a directory
    structure following the top BLAST hit taxonomy ID (taxid) and scientific name.

    For example:

    \b
    taxonomy-sorted-sequences/
    ├── 0-no-hits/
    │   ├── Sample2.fasta
    ├── Foot_and_mouth_disease_virus-12110/
    │   ├── Sample1.fasta
    │   └── Sample2.fasta
    └── Bluetongue_virus-40051/
        ├── Sample1.fasta
        └── Sample3.fasta

    Where the output directory name contains the scientific name of the match with
    all non-word characters replaced with underscore ('_') and the taxid ("{sciname}-{taxid}").

    \033[1mNote:\033[0m Unclassified sequences for each input will be put in the "unclassified-0" directory if there are
    any such sequences. It may be useful to look at these sequences in more detailed with more in-depth searches
    using larger BLAST databases or other tools.
    """

    init_console_logger(verbose)
    blast_tsv_dir = Path(blast_tsv_dir)

    if blast_tsv_sample_name_pattern:
        sample_name_regex = re.compile(blast_tsv_sample_name_pattern)
        logger.debug(f'sample_name_regex={sample_name_regex}')
        sample_blast_tsv: Dict[str, Path] = {sample_name_regex.sub(r'\1', x.name).replace('.tsv', ''): x for x in blast_tsv_dir.glob('*.tsv')}
    else:
        sample_blast_tsv: Dict[str, Path] = {x.name.replace('.tsv', ''): x for x in blast_tsv_dir.glob('*.tsv')}
    logger.info(f'Found {len(sample_blast_tsv)} TSV files in "{blast_tsv_dir}"')
    sample_dfs = blast_tsv_to_df(sample_blast_tsv, top_n_results)
    output_xlsx_report(output_xlsx, sample_dfs, top_n_results)
    if seq_dir:
        fastas, present_fastas = check_fastas(sample_blast_tsv, seq_dir, seq_outdir)
        path_seqoutdir = Path(seq_outdir)
        path_seqoutdir.mkdir(parents=True, exist_ok=True)
        write_seqs_to_taxonomy_dirs(fastas,
                                    path_seqoutdir,
                                    present_fastas,
                                    sample_dfs,
                                    keep_orientation=keep_orientation)


def check_fastas(sample_blast_tsv: Dict[str, Path], seq_dir: str, seq_outdir: str) -> Tuple[Dict[str, Path], Set[str]]:
    seq_dir = Path(seq_dir)
    logger.info(f'FASTA sequence directory provided: "{seq_dir}". Taxonomy sorted sequences will be output'
                f' to "{seq_outdir}".')
    fastas = {x.name.replace('.fasta', ''): x for x in seq_dir.glob('*.fasta')}
    if len(fastas) == 0:
        logger.warning(f'FASTA sequence directory contains no FASTA files matching glob pattern "*.fasta"!')
    # check that all BLAST TSV files have an associated FASTA file
    missing_fastas = set(sample_blast_tsv.keys()) - set(fastas.keys())
    if len(missing_fastas) > 0:
        logger.warning(f'FASTA files missing for the following samples: {missing_fastas}')
    present_fastas = set(sample_blast_tsv.keys()) & set(fastas.keys())
    logger.debug(f'Sorting sequences from following FASTAs: {present_fastas}')
    return fastas, present_fastas


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
