# coding: utf-8
import logging
import re
from pathlib import Path
from typing import Tuple, Iterator, List, Union, Dict, Set, Mapping

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from blast2xl.util import invert_dict

logger = logging.getLogger(__name__)

BLAST_COLUMNS = [('qaccver', 'Query', 'category'),
                 ('saccver', 'Subject', 'category'),
                 ('pident', 'Percent_Identity', float),
                 ('length', 'Alignment_Length', 'uint32'),
                 ('mismatch', '#_Mismatches', 'uint32'),
                 ('gapopen', '#_Gap_Openings', 'uint32'),
                 ('qstart', 'Query_Start', 'uint32'),
                 ('qend', 'Query_End', 'uint32'),
                 ('sstart', 'Subject_Start', 'uint32'),
                 ('send', 'Subject_End', 'uint32'),
                 ('evalue', 'e_Value', float),
                 ('bitscore', 'Bitscore', float),
                 ('qlen', 'Query_Length', 'uint32'),
                 ('slen', 'Subject_Length', 'uint32'),
                 # ('sstrand', 'Subject_Strand', 'category'),
                 ('stitle', 'Subject_Title', str),
                 ('staxid', 'Subject_taxid', 'uint32'),
                 ('ssciname', 'Subject_Sciname', 'category')]


def read_blast_tsv(blast_tsv_path: Union[str, Path]) -> pd.DataFrame:
    return pd.read_table(blast_tsv_path,
                         header=None,
                         names=[y for x, y, z in BLAST_COLUMNS],
                         dtype={y: z for x, y, z in BLAST_COLUMNS})


def get_col_widths(df: pd.DataFrame, index: bool = False, pad_width: int = 2) -> Iterator[int]:
    """Calculate column widths based on column headers and contents

    Supports multi-index
    """
    if index:
        for i, idx_name in enumerate(df.index.names):
            idx_max = max([len(str(s[i] if len(s) > 1 else s)) for s in df.index.values] + [len(str(idx_name))])
            yield idx_max + pad_width
    for c in df.columns:
        # get max length of column contents and length of column header
        yield np.max([df[c].astype(str).str.len().max() + 1, len(c) + 1]) + pad_width


def write_excel(name_dfs: List[Tuple[str, pd.DataFrame]],
                output_dest: str,
                output_df_index: bool = False,
                sheet_name_index: bool = True,
                freeze_panes: Tuple[int, int] = (1, 1)) -> None:
    if not output_dest.endswith('.xlsx'):
        output_dest += '.xlsx'
    logger.info('Starting to write Pandas DataFrames to worksheets in XLSX workbook ("{output_dest}")')
    with pd.ExcelWriter(output_dest, engine='xlsxwriter') as writer:
        forbidden_characters = re.compile(r'[\\:/?*\[\]]+')
        idx = 1
        for name_df in name_dfs:
            if not isinstance(name_df, (list, tuple)):
                logger.error(f'Input "{name_df}" is not a list or tuple (type="{type(name_df)}"). Skipping...')
                continue
            sheetname, df = name_df
            fixed_sheetname = forbidden_characters.sub('_', sheetname)
            # fixed max number of characters in sheet name due to compatibility
            if sheet_name_index:
                max_chars = 28
                fixed_sheetname = f'{idx}_{fixed_sheetname[:max_chars]}'
            else:
                max_chars = 31
                fixed_sheetname = fixed_sheetname[:max_chars]

            if len(fixed_sheetname) > max_chars:
                logger.warning(f'Sheetname "{fixed_sheetname}" is >= {max_chars} characters so may be truncated '
                               f'(n={len(fixed_sheetname)})')

            logger.info(f'Writing table to Excel sheet "{fixed_sheetname}"')
            df.to_excel(writer,
                        sheet_name=fixed_sheetname,
                        index=output_df_index,
                        freeze_panes=freeze_panes)
            worksheet = writer.book.get_worksheet_by_name(fixed_sheetname)
            for i, width in enumerate(get_col_widths(df, index=output_df_index)):
                worksheet.set_column(i, i, width)
            idx += 1
    logger.info('Done writing worksheets to spreadsheet "%s".', output_dest)


def write_seqs_to_taxonomy_dirs(fastas: Dict[str, Path],
                                path_seqoutdir: Path,
                                present_fastas: Set[str],
                                sample_dfs: Dict[str, pd.DataFrame],
                                keep_orientation: bool = False):
    for sample in present_fastas:
        df = sample_dfs[sample]
        df_top: pd.DataFrame = df.groupby('Query').first()
        qid_dirname_series: pd.Series = df_top.Subject_Sciname.str.replace(r'\W', '_') \
                                        + '-' + df_top.Subject_taxid.astype(str)
        query_dirname: Mapping[str, str] = qid_dirname_series.to_dict()
        dirname_queries = invert_dict(query_dirname)
        recs: Dict[str, SeqRecord] = {r.id: r for r in SeqIO.parse(str(fastas[sample]), 'fasta')}
        missing_recs = set(recs.keys()) - set(query_dirname.keys())
        if not keep_orientation:
            need_revcomp = set(df_top.index.values[df_top.Subject_Start > df_top.Subject_End])
            out = {}
            for rid, r in recs.items():
                if rid in need_revcomp:
                    r_desc = r.description
                    r_name = r.name
                    r = r.reverse_complement()
                    r.id = rid
                    r.name = r_name
                    r.description = f'{r_desc}|REVCOMP'
                if rid not in missing_recs:
                    rid_blast = df_top.loc[rid, :].to_dict()
                    r.description += (f'|TOP_ACC="{rid_blast["Subject"]}"'
                                      f'|TOP_%ID={rid_blast["Percent_Identity"]:.2f}'
                                      f'|TOP_ALN_LEN={rid_blast["Alignment_Length"]}'
                                      f'|TOP_NAME="{rid_blast["Subject_Title"]}"')
                    # r.description += f'|{rid_blast}'
                out[rid] = r
            recs = out

        if len(missing_recs) > 0:
            unclassified_dir = path_seqoutdir / sample / '0-no-hits'
            unclassified_dir.mkdir(parents=True, exist_ok=True)
            sample_seqout = unclassified_dir / f'{sample}.fasta'
            n_written = SeqIO.write([recs[r] for r in missing_recs], sample_seqout, 'fasta')
            logger.info(f'Wrote {n_written} sequences to "{sample_seqout}"')
        for dirname, queries in dirname_queries.items():
            path_sciname = path_seqoutdir / sample / dirname
            path_sciname.mkdir(parents=True, exist_ok=True)
            sample_seqout = path_sciname / f'{sample}.fasta'
            n_written = SeqIO.write([recs[q] for q in queries], sample_seqout, 'fasta')
            logger.info(f'Wrote {n_written} sequences to "{sample_seqout}"')
