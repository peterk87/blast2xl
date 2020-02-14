# coding: utf-8
"""Main module."""
import logging
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

from blast2xl.io import read_blast_tsv, write_excel

logger = logging.getLogger(__name__)


def blast_tsv_to_df(sample_blast_tsv: Dict[str, Path], top_n_results: int = -1) -> Dict[str, pd.DataFrame]:
    sample_dfs = {}
    for sample, tsv_path in sample_blast_tsv.items():
        logger.debug(f'Parsing sample "{sample}" tabular BLAST result into DataFrame ({tsv_path})')
        df = read_blast_tsv(tsv_path)
        logger.debug(f'Parsed sample "{sample}" tabular BLAST result into DataFrame with {df.shape[0]} rows')
        df['Sample'] = sample
        df.sort_values(['Query', 'Bitscore'], ascending=[True, False], inplace=True)
        sample_dfs[sample] = df.groupby('Query').head(top_n_results) if top_n_results > 0 else df
    return sample_dfs


def output_xlsx_report(output_path: str,
                       sample_blast: Dict[str, pd.DataFrame],
                       top_n_results: int = -1):
    df_concat: pd.DataFrame = pd.concat((df for sample, df in sample_blast.items()),
                                        sort=False,
                                        ignore_index=True) \
        .sort_values(['Sample', 'Query'], ascending=True)
    subj_accessions = df_concat['Subject'].astype(str)
    if subj_accessions.str.match(r'^[A-Z]{1,2}\d{5,}(\.\d+)?$').all():
        df_concat['Subject_URL'] = 'https://www.ncbi.nlm.nih.gov/nuccore/' + subj_accessions
    taxids = df_concat['Subject_taxid'].astype(str)
    df_concat['Subject_taxid_URL'] = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=' + taxids
    df_top_results_per_query_seq = df_concat \
        .groupby(['Sample', 'Query']) \
        .first() \
        .reset_index() \
        .sort_values(['Sample', 'Query'], ascending=True)
    all_results_sheetname = 'All BLAST Results' if top_n_results <= 0 else f'Top {top_n_results} BLAST Results'

    summaries = []
    for sample, df in sample_blast.items():
        logger.debug(f'{sample}: {df.shape}')
        df_copy = df.copy()
        df_copy.Subject = df_copy.Subject.astype(str)
        df_summary = df_copy.reset_index().groupby('Subject').apply(subj_cov).head(5)
        df_summary.sort_values(['total_aligned_positions', 'avg_pid'], ascending=False, inplace=True)
        df_summary.avg_pid /= 10000
        df_summary['p_aligned'] = df_summary.total_aligned_positions / df_summary.subject_length
        df_summary['sample'] = sample
        summaries.append(df_summary.reset_index().set_index('sample'))

    df_summary = pd.concat(summaries)
    write_excel([('Summary', df_summary),
                 ('Top BLAST Results', df_top_results_per_query_seq.set_index(['Sample', 'Query'])),
                 (all_results_sheetname, df_concat.set_index(['Sample', 'Query']))],
                output_dest=output_path,
                output_df_index=True,
                sheet_name_index=False,
                freeze_panes=(1, 2))


def subj_cov(dfg: pd.DataFrame) -> pd.Series:
    """Summarize subject sequence alignment info

    Total alignment coverage of query sequences for a given subject, average
    %ID and other stats to determine the top overall subject sequence.
    """
    logger.debug(f'dfg subject "{dfg.index.unique()}"| shape={dfg.shape}')
    if dfg.shape[0] == 0:
        return pd.Series()
    arrlen = dfg.Subject_Length.values[0]
    # forward and reverse subject coverage
    cov_mask = np.zeros((2, arrlen), dtype=bool)

    starts = dfg.Subject_Start.values
    ends = dfg.Subject_End.values
    for start, end in zip(starts, ends):
        direction = 0
        if start > end:
            direction = 1
            start, end = end, start
        cov_mask[direction, start:(end + 1)] = True
    forward_cov = cov_mask[0].sum()
    reverse_cov = cov_mask[1].sum()
    total_cov = (cov_mask[0] | cov_mask[1]).sum()
    fwd_rev_cov_overlap = (cov_mask[0] & cov_mask[1]).sum()
    N = dfg.index.size
    pids = dfg.Percent_Identity.values
    alns = dfg.Alignment_Length.values
    pid_aln = int(np.sum(pids * alns) / np.sum(alns) * 100)

    return pd.Series([total_cov,
                      arrlen,
                      arrlen - total_cov,
                      forward_cov,
                      reverse_cov,
                      fwd_rev_cov_overlap,
                      N,
                      pid_aln],
                     index=['total_aligned_positions',
                            'subject_length',
                            'unaligned_positions',
                            'forward_aligned_positions',
                            'reverse_aligned_positions',
                            'forward_reverse_aligned_positions_overlap',
                            'n_alignments',
                            'avg_pid'])
