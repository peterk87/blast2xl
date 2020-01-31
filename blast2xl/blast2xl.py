# coding: utf-8
"""Main module."""
import logging
from pathlib import Path
from typing import Dict

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
    write_excel([('Top BLAST Results', df_top_results_per_query_seq.set_index(['Sample', 'Query'])),
                 (all_results_sheetname, df_concat.set_index(['Sample', 'Query']))],
                output_dest=output_path,
                output_df_index=True,
                sheet_name_index=False,
                freeze_panes=(1, 2))
