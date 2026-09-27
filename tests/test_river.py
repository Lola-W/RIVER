"""Scientific contract checks using small synthetic alignments.

Run through the supplied scheduler example on a cluster. Real study BAMs are
deliberately not distributed with these fixtures.
"""
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import pysam

from river.inputs import load_samples, load_variants
from river.pileup import parse_read_bases, wilson_interval
from river.selection import passes_read_qc, canonical_umi_pair, classify_f1r2_f2r1


def write_tsv(path, fields, rows):
    with path.open('w', newline='') as handle:
        out = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        out.writeheader()
        out.writerows(rows)


def read_tsv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def record(name='read_AC+GT', flag=99, start=10, chrom=0,
           alt=False, tags=(), cigar='100M'):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.flag = flag
    read.reference_id = chrom
    read.reference_start = start
    read.mapping_quality = 60
    read.cigarstring = cigar
    read.next_reference_id = chrom
    read.next_reference_start = 250
    read.template_length = 340
    query = list('A' * 100)
    if alt:
        query[49 - start] = 'G'
    read.query_sequence = ''.join(query)
    read.query_qualities = pysam.qualitystring_to_array('I' * 100)
    read.set_tags(list(tags))
    return read


class NumericalTests(unittest.TestCase):
    def test_snv_counts_skip_long_indels_and_special_symbols(self):
        text = '^].+123' + 'A' * 123 + ',-12' + 'T' * 12 + 'AaTtNn*#><$'
        obs = parse_read_bases(text, 'A')
        self.assertEqual([(b, rev) for b, rev in obs if b in ('A', 'C', 'G', 'T')],
                         [('A', False), ('A', True), ('A', False), ('A', True),
                          ('T', False), ('T', True)])
        self.assertEqual(len(obs), 12)

    def test_malformed_pileup_fails(self):
        for text in ('^', '.+3AA', '.+A', '.!', '.-100AC'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_read_bases(text, 'A')

    def test_wilson_matches_reference_values(self):
        vaf, low, high = wilson_interval(20, 55)
        self.assertAlmostEqual(vaf, 20 / 55, places=14)
        self.assertAlmostEqual(low, 0.24930531530307248, places=14)
        self.assertAlmostEqual(high, 0.49577238500158743, places=14)

    def test_threshold_boundary_does_not_pass(self):
        _, low, _ = wilson_interval(2, 549)
        self.assertAlmostEqual(low, 0.00099960765059306, places=15)
        self.assertFalse(low > 0.001)

    def test_zero_depth_is_undefined(self):
        for value in wilson_interval(0, 0):
            self.assertTrue(value is None or math.isnan(value))


class SelectionRuleTests(unittest.TestCase):
    def test_cigar_length_four_is_literal(self):
        for cigar, expected in [('4S96M', False), ('5S95M', True),
                                ('4M1I95M', False), ('100M', True), ('96M4S', True)]:
            with self.subTest(cigar=cigar):
                self.assertEqual(passes_read_qc(record(cigar=cigar)), expected)

    def test_qc_rules_and_duplicate_family_eligibility(self):
        self.assertTrue(passes_read_qc(record(flag=99 | 1024)))
        for flag in (99 | 256, 99 | 512, 99 | 2048, 97):
            self.assertFalse(passes_read_qc(record(flag=flag)))
        self.assertFalse(passes_read_qc(record(tags=[('DT', 'SQ')])))

    def test_complementary_umi_order(self):
        first = record(name='first_AC+GT', flag=99)
        second = record(name='second_GT+AC', flag=83)
        self.assertEqual(canonical_umi_pair(first), canonical_umi_pair(second))
        self.assertEqual(classify_f1r2_f2r1(first), 1)
        self.assertEqual(classify_f1r2_f2r1(second), 2)


class InputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_source_only_and_repeated_origins(self):
        sample_file = self.root / 'samples.tsv'
        # Input parsing is independent of opening/indexing the alignment.
        (self.root / 'reads.bam').touch()
        write_tsv(sample_file, ['sample_id', 'method', 'type', 'bam'], [
            dict(sample_id='source', method='Caller', type='bulk', bam='.'),
            dict(sample_id='query', method='IlluminaSeq', type='bulk', bam='reads.bam')])
        samples = load_samples(sample_file)
        self.assertIsNone(samples['source'].bam)
        self.assertEqual(samples['query'].bam, (self.root / 'reads.bam').resolve())
        variant_file = self.root / 'variants.tsv'
        row = dict(chrom='chr1', pos=50, ref='A', alt='G', source_sample_id='source')
        write_tsv(variant_file, list(row), [row, row])
        self.assertEqual(len(load_variants(variant_file, samples)), 1)

    def test_inconsistent_method_type_rejected(self):
        path = self.root / 'samples.tsv'
        write_tsv(path, ['sample_id', 'method', 'type', 'bam'], [
            dict(sample_id='a', method='Custom', type='bulk', bam='.'),
            dict(sample_id='b', method='Custom', type='single_cell', bam='.')])
        with self.assertRaises(ValueError):
            load_samples(path)

    def test_conflicting_preset_type_rejected(self):
        path = self.root / 'samples.tsv'
        write_tsv(path, ['sample_id', 'method', 'type', 'bam'], [
            dict(sample_id='a', method='UDSeq', type='single_cell', bam='.')])
        with self.assertRaises(ValueError):
            load_samples(path)


class IntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.samtools = os.environ.get('RIVER_TEST_SAMTOOLS') or shutil.which('samtools')
        if not cls.samtools:
            raise unittest.SkipTest('Set RIVER_TEST_SAMTOOLS to run BAM integration tests')

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.reference = self.root / 'reference.fa'
        self.reference.write_text('>chr1\n' + 'A' * 500 + '\n>chr2\n' + 'A' * 500 + '\n')
        pysam.faidx(str(self.reference))

    def bam(self, name, records):
        path = self.root / (name + '.bam')
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': 'chr1', 'LN': 500}, {'SN': 'chr2', 'LN': 500}]}
        with pysam.AlignmentFile(str(path), 'wb', header=header) as out:
            for read in sorted(records, key=lambda r: (r.reference_id, r.reference_start)):
                out.write(read)
        pysam.index(str(path))
        return path.name

    def run_river(self, samples, variants, name='results', expect_success=True):
        sample_path = self.root / (name + '.samples.tsv')
        fields = ['sample_id', 'method', 'type', 'bam']
        fields += sorted(set().union(*(row.keys() for row in samples)) - set(fields))
        write_tsv(sample_path, fields, samples)
        variant_path = self.root / (name + '.variants.tsv')
        write_tsv(variant_path, ['chrom', 'pos', 'ref', 'alt', 'source_sample_id', 'callset'], variants)
        out = self.root / name
        proc = subprocess.run([sys.executable, '-m', 'river', 'run', '--variants', str(variant_path),
                               '--samples', str(sample_path), '--reference', str(self.reference),
                               '--samtools', self.samtools, '--outdir', str(out)],
                              capture_output=True, text=True)
        if expect_success:
            self.assertEqual(proc.returncode, 0, proc.stdout + '\n' + proc.stderr)
            self.assertEqual(json.loads((out / 'run.json').read_text())['status'], 'complete')
        else:
            self.assertNotEqual(proc.returncode, 0, proc.stdout + '\n' + proc.stderr)
        return out, proc

    @staticmethod
    def candidate(source='own', pos=50, ref='A', alt='G'):
        return dict(chrom='chr1', pos=pos, ref=ref, alt=alt, source_sample_id=source, callset='test')

    def test_method_aggregation_and_denominators(self):
        definitions = [
            ('own', 'IlluminaSeq', 'bulk', 10, 0),
            ('replicate', 'IlluminaSeq', 'bulk', 100, 20),
            ('duplex', 'UDSeq', 'duplex', 100, 5),
            ('cell1', 'PTA-Seq', 'single_cell', 1, 1),
            ('cell2', 'PTA-Seq', 'single_cell', 1, 1)]
        samples = []
        for name, method, kind, n, alt in definitions:
            bam = self.bam(name, [record(name=f'{name}_{i}', flag=0, alt=i < alt) for i in range(n)])
            samples.append(dict(sample_id=name, method=method, type=kind, bam=bam,
                                bam_state='selected' if kind == 'duplex' else ''))
        samples.append(dict(sample_id='source_only', method='OtherCaller', type='bulk', bam='.'))
        variants = [self.candidate(), self.candidate(), self.candidate(pos=400),
                    self.candidate(source='cell1'), self.candidate(source='source_only')]
        out, _ = self.run_river(samples, variants)
        rows = read_tsv(out / 'read_support.tsv')
        self.assertEqual(len(rows), 10)
        evidence = {(r['sample_id'], int(r['pos'])): r for r in rows}
        self.assertEqual(int(evidence[('replicate', 50)]['alt_count']), 20)
        self.assertEqual(int(evidence[('duplex', 50)]['ref_count']), 95)
        for sample in definitions:
            self.assertEqual(int(evidence[(sample[0], 400)]['depth']), 0)
        results = {(r['source_sample_id'], int(r['pos'])): r for r in read_tsv(out / 'variant_support.tsv')}
        self.assertEqual(len(results), 4)
        self.assertEqual(int(results[('own', 50)]['n_other_supporting_methods']), 2)
        self.assertEqual(int(results[('cell1', 50)]['n_other_supporting_methods']), 2)
        self.assertEqual(int(results[('source_only', 50)]['n_other_supporting_methods']), 3)
        self.assertEqual(int(results[('own', 400)]['n_other_supporting_methods']), 0)
        self.assertEqual(json.loads(results[('own', 50)]['supporting_cell_counts_by_method']),
                         {'PTA-Seq': 2})
        self.assertEqual(json.loads(results[('own', 50)]['supporting_cell_classes_by_method']),
                         {'PTA-Seq': '>=2'})
        summary = {r['source_sample_id']: r for r in read_tsv(out / 'summary.tsv')}
        self.assertEqual(int(summary['own']['n_variants']), 2)
        self.assertEqual(int(summary['own']['n_support_0']), 1)
        self.assertEqual(int(summary['own']['n_support_ge2']), 1)
        # Reusing an output directory must fail without replacing a completed run.
        saved = (out / 'summary.tsv').read_bytes()
        self.run_river(samples, variants, expect_success=False)
        self.assertEqual((out / 'summary.tsv').read_bytes(), saved)

    def test_raw_presets_and_custom_selection(self):
        nano = [record('a_AC+GT', flag=99, alt=True), record('b_GT+AC', flag=83),
                record('one_AC+CC', flag=99), record('bad1_AA+AA', flag=99, tags=[('DT', 'SQ')]),
                record('bad2_AA+AA', flag=83)]
        ppm = [record('pass', flag=0, alt=True, tags=[('st', 'MIXED'), ('et', 'MIXED')]),
               record('fail', flag=0, tags=[('st', 'MIXED'), ('et', 'PLUS')]),
               record('dup', flag=1024, tags=[('st', 'MIXED'), ('et', 'MIXED')])]
        hidef = [record('movie/3/ccs/fwd', flag=0, alt=True),
                 record('movie/3/ccs/rev', flag=16, chrom=1),
                 record('movie/3/ccs/other', flag=0), record('movie/4/ccs/fwd', flag=0)]
        custom = [record('c1', flag=0, alt=True, tags=[('MI', 'one'), ('DS', 'F')]),
                  record('c2', flag=16, tags=[('MI', 'one'), ('DS', 'R')]),
                  record('c3', flag=0, tags=[('MI', 'two'), ('DS', 'F')]),
                  record('c4', flag=16, chrom=1, tags=[('MI', 'two'), ('DS', 'R')])]
        samples = []
        for method, reads in [('UDSeq', nano), ('NanoSeq', nano), ('ppmSeq', ppm), ('HiDEF-seq', hidef)]:
            samples.append(dict(sample_id=method, method=method, type='duplex',
                                bam=self.bam(method, reads)))
        samples.append(dict(sample_id='custom', method='MyDuplex', type='duplex',
                            bam=self.bam('custom', custom), strand_source='tag:DS',
                            fwd_pattern='^F$', rev_pattern='^R$', family_tag='MI'))
        out, _ = self.run_river(samples, [self.candidate(source='UDSeq')])
        evidence = {r['sample_id']: r for r in read_tsv(out / 'read_support.tsv')}
        expected = {'UDSeq': (1, 1), 'NanoSeq': (1, 1), 'ppmSeq': (0, 1),
                    'HiDEF-seq': (1, 1), 'custom': (1, 1)}
        for name, (ref, alt) in expected.items():
            with self.subTest(method=name):
                self.assertEqual(int(evidence[name]['ref_count']), ref)
                self.assertEqual(int(evidence[name]['alt_count']), alt)

    def test_bad_reference_and_missing_index_are_errors(self):
        name = self.bam('input', [record(flag=0)])
        samples = [dict(sample_id='own', method='IlluminaSeq', type='bulk', bam=name)]
        self.run_river(samples, [self.candidate(ref='C')], name='wrong_ref', expect_success=False)
        (self.root / (name + '.bai')).unlink()
        self.run_river(samples, [self.candidate()], name='missing_index', expect_success=False)


    def test_orientation_labels_do_not_gate_bulk_counts(self):
        reads = [record('f', flag=99, alt=True), record('r', flag=83, alt=True),
                 record('unpaired', flag=0, alt=True), record('ref', flag=147)]
        samples = [dict(sample_id='own', method='IlluminaSeq', type='bulk',
                        bam=self.bam('orientations', reads), strand_source='sam_flags',
                        fwd_pattern='F1R2', rev_pattern='F2R1')]
        out, _ = self.run_river(samples, [self.candidate()])
        row = read_tsv(out / 'read_support.tsv')[0]
        self.assertEqual([int(row[key]) for key in ('alt_forward', 'alt_reverse', 'alt_unclassified')], [1, 1, 1])
        self.assertEqual(int(row['ref_forward']), 1)
        for allele in ('ref', 'alt'):
            self.assertEqual(int(row[allele + '_count']),
                             sum(int(row[allele + '_' + label]) for label in ('forward', 'reverse', 'unclassified')))

    def test_tag_orientation_missing_and_ambiguous_labels(self):
        reads = [record('f', flag=0, alt=True, tags=[('DS', 'F')]),
                 record('r', flag=0, alt=True, tags=[('DS', 'R')]), record('unknown', flag=0, alt=True)]
        samples = [dict(sample_id='own', method='CustomBulk', type='bulk',
                        bam=self.bam('tags', reads), strand_source='tag:DS', fwd_pattern='F', rev_pattern='R')]
        out, _ = self.run_river(samples, [self.candidate()])
        row = read_tsv(out / 'read_support.tsv')[0]
        self.assertEqual([int(row[key]) for key in ('alt_count', 'alt_forward', 'alt_reverse', 'alt_unclassified')], [3, 1, 1, 1])
        samples[0]['bam'] = self.bam('ambiguous', [record('both', flag=0, tags=[('DS', 'FR')])])
        self.run_river(samples, [self.candidate()], name='ambiguous', expect_success=False)

    def test_multiple_alts_and_first_reference_position(self):
        third = record('third', flag=0)
        third.query_sequence = 'A' * 39 + 'T' + 'A' * 60
        third.query_qualities = pysam.qualitystring_to_array('I' * 100)
        edge = record('edge', flag=0, start=0, cigar='1M')
        edge.query_sequence = 'A'
        edge.query_qualities = pysam.qualitystring_to_array('I')
        samples = [dict(sample_id='own', method='IlluminaSeq', type='bulk',
                        bam=self.bam('alleles', [record('ref', flag=0), record('alt', flag=0, alt=True), third, edge]))]
        out, _ = self.run_river(samples, [self.candidate(), self.candidate(alt='T'), self.candidate(pos=1)])
        rows = {(int(r['pos']), r['alt']): r for r in read_tsv(out / 'read_support.tsv')}
        for alt in ('G', 'T'):
            self.assertEqual(int(rows[(50, alt)]['ref_count']), 1)
            self.assertEqual(int(rows[(50, alt)]['alt_count']), 1)
            self.assertEqual(int(rows[(50, alt)]['depth']), 2)
        self.assertEqual(int(rows[(1, 'G')]['ref_count']), 1)

    def test_same_method_counts_are_not_pooled(self):
        samples = [dict(sample_id='own', method='Caller', type='bulk', bam='.')]
        for dataset in ('a', 'b'):
            reads = [record(f'{dataset}_{i}', flag=0, alt=i == 0) for i in range(200)]
            samples.append(dict(sample_id=dataset, method='RepeatedBulk', type='bulk', bam=self.bam(dataset, reads)))
        out, _ = self.run_river(samples, [self.candidate()])
        self.assertTrue(all(int(r['qualifying_support']) == 0 for r in read_tsv(out / 'read_support.tsv')))
        self.assertEqual(int(read_tsv(out / 'variant_support.tsv')[0]['n_other_supporting_methods']), 0)
        self.assertGreater(wilson_interval(2, 400)[1], 0.001)

    def test_raw_scheme_mismatch_errors_but_valid_zero_is_allowed(self):
        specs = [
            ('ppmSeq', {}, [record('missing_tags', flag=0)],
             [record('unmixed', flag=0, tags=[('st', 'PLUS'), ('et', 'PLUS')])]),
            ('HiDEF-seq', {}, [record('unrecognized_name', flag=0)],
             [record('movie/3/ccs/fwd', flag=0)]),
            ('NewDuplex', dict(strand_source='tag:DS', fwd_pattern='^F$', rev_pattern='^R$', family_tag='MI'),
             [record('unknown', flag=0, tags=[('MI', 'one'), ('DS', 'X')])],
             [record('one_sided', flag=0, tags=[('MI', 'one'), ('DS', 'F')])])]
        for index, (method, options, invalid, valid_zero) in enumerate(specs):
            with self.subTest(method=method):
                row = dict(sample_id='own', method=method, type='duplex',
                           bam=self.bam(f'bad_scheme_{index}', invalid), **options)
                out, _ = self.run_river([row], [self.candidate()], name=f'bad_scheme_{index}', expect_success=False)
                self.assertEqual(json.loads((out / 'run.json').read_text())['status'], 'failed')
                row['bam'] = self.bam(f'valid_zero_{index}', valid_zero)
                out, _ = self.run_river([row], [self.candidate()], name=f'valid_zero_{index}')
                self.assertEqual(int(read_tsv(out / 'read_support.tsv')[0]['depth']), 0)
                row['bam'] = self.bam(f'empty_{index}', [])
                out, _ = self.run_river([row], [self.candidate()], name=f'empty_{index}')
                self.assertEqual(int(read_tsv(out / 'read_support.tsv')[0]['depth']), 0)

    def test_failed_mpileup_cannot_publish_zero_results(self):
        real_tool = self.samtools
        fake = self.root / 'failed_samtools'
        fake.write_text('#!/usr/bin/env bash\nif [[ "$1" == "--version" ]]; then echo test; exit 0; fi\necho intentional_failure >&2\nexit 9\n')
        fake.chmod(0o755)
        self.samtools = str(fake)
        try:
            samples = [dict(sample_id='own', method='IlluminaSeq', type='bulk',
                            bam=self.bam('valid', [record(flag=0)]))]
            out, _ = self.run_river(samples, [self.candidate()], expect_success=False)
            self.assertEqual(json.loads((out / 'run.json').read_text())['status'], 'failed')
            self.assertFalse((out / 'read_support.tsv').exists())
            self.assertFalse((out / 'summary.tsv').exists())
        finally:
            self.samtools = real_tool


if __name__ == '__main__':
    unittest.main()
