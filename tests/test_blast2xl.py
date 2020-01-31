#!/usr/bin/env python

"""Tests for `blast2xl` package."""

from os.path import abspath
from pathlib import Path

from click.testing import CliRunner

from blast2xl import cli


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert 'blast2xl: BLAST XLSX Report Creator' in help_result.stdout

    result = runner.invoke(cli.main)
    assert result.exit_code == 2
    assert 'Missing option' in result.output

    blast_tsv_dir = abspath('tests/data/blast_tsv')
    fasta_dir = abspath('tests/data/fastas')
    with runner.isolated_filesystem():
        excel_report = 'blast-report.xlsx'
        seq_outdir = 'seq-outdir'
        result = runner.invoke(cli.main, ['--blast-tsv-dir', blast_tsv_dir,
                                          '--blast-tsv-sample-name-pattern', r'^blastn-(.+)-vs-nt.*',
                                          '--seq-dir', fasta_dir,
                                          '--top-n-results', 5,
                                          '-o', excel_report,
                                          '-O', seq_outdir,
                                          '-vvv'])

        assert result.exit_code == 0
        path_seq_outdir = Path(seq_outdir)
        assert path_seq_outdir.exists()
        output_fastas = list(path_seq_outdir.glob('**/*.fasta'))
        assert len(output_fastas) > 2

        fasta_path = path_seq_outdir / 'FMDV' / 'Foot_and_mouth_disease_virus___type_O-12118' / 'FMDV.fasta'
        assert fasta_path.exists()
        assert fasta_path.stat().st_size > 0
        assert Path(excel_report).exists()
        assert Path(excel_report).stat().st_size > 0

