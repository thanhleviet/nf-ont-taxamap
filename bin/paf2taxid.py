#!/usr/bin/env python

import click
from loguru import logger
import gzip
from ete3 import NCBITaxa
import pandas as pd


class Alignment(object):
    """
    Create alignment object based on PAF line.
    """
    def __init__(self, line):
        self.line = line
        self.query_name = 0
        self.query_length = 0
        self.query_start = 0
        self.query_end = 0
        self.ref_name = None
        self.ref_length = 0
        self.ref_start = 0
        self.ref_end = 0
        self.matched_mapping_bases = 0.0
        self.total_mapping_bases = 0.0
        self.score = 0.0
        self.taxid = None
        self.ref_accession = None
        self._parse_line()

    def _parse_line(self):
        cols = self.line.decode().strip().split("\t")
        if len(cols) >= 12:
            self.query_name, self.query_length, self.query_start, self.query_end, _, self.ref_name, self.ref_length, self.ref_start, self.ref_end, self.matched_mapping_bases, self.total_mapping_bases = cols[:11]
            self.aligned_query_length = int(self.query_end) - int(self.query_start)
            self.ref_name = self.ref_name.replace("kraken:taxid|", "")
            
            self.aligned_coverage = int(self.matched_mapping_bases) / int(self.query_length)* 100
            self.aligned_coverage = round(self.aligned_coverage, 2)

            self.aligned_identity = float(self.matched_mapping_bases) / float(self.total_mapping_bases) * 100
            self.aligned_identity = round(self.aligned_identity, 2)
            
            self.read_coverage = float(self.aligned_query_length) / float(self.query_length) * 100
            self.read_coverage = round(self.read_coverage, 2)

            self.score = 2*((self.aligned_coverage * self.aligned_identity)/(self.aligned_coverage + self.aligned_identity))
            self.score = round(self.score, 2)

            self.taxid, self.ref_accession = self.ref_name.split("|")
        else:
            logger.error(f"Incorrect line: {cols}")


def write_taxonomy_read(file_name, collection):
    try:
        with open(file_name, 'w') as _fh:
            _fh.write("read,taxid,read_length,closet_reference,read_coverage,alignment_coverage,alignment_identity,score\n")
            for k, v in collection.items():
                _fh.write(
                    f"{k},{collection[k]['taxid']},{collection[k]['read_length']},{collection[k]['ref_accession']},{collection[k]['read_coverage']},{collection[k]['aligned_coverage']},{collection[k]['identity']},{collection[k]['score']}\n")
    except IOError as e:
        logger.error("Could not write")
        logger.error(e)


def write_taxonomy_report(file_name, collection):
    try:
        with open(file_name, 'w') as _fh:
            _fh.write(
                "taxid,count,kingdom,phylum,class,order,family,genus,species,strain\n")
            for line in collection:
                _fh.write(f"{line}\n")
    except IOError as e:
        logger.error("Could not write")
        logger.error(e)


def write_metaphlan_like_report(file_name, report_file):
    try:
        ncbi = NCBITaxa()
        report = pd.read_csv(report_file)
        taxonomy = [
            "kingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
            "strain"
        ]

        _list = []
        _output = []
        for rank in taxonomy:
            _list.append(rank)
            _agg = report.groupby(_list)['count'].agg('sum').reset_index()
            df_to_list = _agg.values.tolist()
            _output.append(df_to_list)
        with open(file_name, 'w') as _fh:
            _fh.write("#mpa_v3\n")
            _fh.write("#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n")
            for x in _output:
                for k in x:
                    k_ = "|".join([m for m in k if isinstance(m, str)])
                    i_ = [m for m in k if isinstance(m, int)][0]
                    splitted_taxon = k_.rsplit("|", 1)
                    i = 0
                    if len(splitted_taxon) > 1:
                        i = 1
                    _taxid_name = [splitted_taxon[i].split("__")[1].replace('~~', ',')]
                    taxid = list(ncbi.get_name_translator(_taxid_name).values())
                    _line = (
                        f"{k_}\t{taxid[0][0]}\t{i_}\t{''}\n")
                    _fh.write(_line)

    except IOError as e:
        logger.error("Could not write")
        logger.error(e)


def swap_rank_dict(rank_dict: dict()):
    _rank = {}
    for k, v in rank_dict.items():
        _rank[v] = k
    return _rank


def print_taxid(taxid):
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    lineage_name = ncbi.get_taxid_translator(lineage)
    rank = ncbi.get_rank(lineage)
    rank = swap_rank_dict(rank)

    taxonomy_rank = ["superkingdom",
                    "kingdom",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species",
                    "strain"
                    ]

    _taxonomy = {}
    for taxon in taxonomy_rank:
        if taxon in rank.keys():
            _value = lineage_name[rank[taxon]]
            if taxon == "kingdom" and _value != "Fungi":
                _value = lineage_name[rank["superkingdom"]]
        if taxon != "superkingdom":
            _taxonomy[taxon] = _value
    return _taxonomy


def print_metaphlan_like_report(taxid):
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(taxid)
    lineage_name = ncbi.get_taxid_translator(lineage)
    rank = ncbi.get_rank(lineage)
    rank = swap_rank_dict(rank)

    taxonomy_rank_dict = {
        "superkingdom": "k__",
        "kingdom": "k__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
        "strain": "t__"
    }
    _taxonomy = {}

    for taxon in taxonomy_rank_dict.keys():
        if taxon in rank.keys():
            _value = f"{taxonomy_rank_dict[taxon]}{lineage_name[rank[taxon]].replace(',','~~')}"
            # logger.info(f"{taxon} - {lineage_name[rank[taxon]]}")
            # logger.info(f"{taxon}-{_value}")
            if taxon == "kingdom" and _value != "Fungi":
                _value = f"{taxonomy_rank_dict['superkingdom']}{lineage_name[rank['superkingdom']]}"
        if taxon != "superkingdom":
            _taxonomy[taxon] = _value
            _value = ""
    return _taxonomy


@click.command()
@click.argument('paf_file')
@click.option('--output', '-o', help="Output name", required=True)
@click.option('--min_read_length', '-l', help="Min read length", default=100)
@click.option('--score', '-c', help="Score threshold to exclude false positive classification [score = 2*(coverage*identity)/(coverage+identity)*100]", default=35.0)
def main(min_read_length, score, paf_file, output=None):
    """
    Parsing a paf file into a readable taxonomy ID format \n
    thanh.le-viet@quadram.ac.uk
    """
    logger.info(f"Processing file: {paf_file}")
    try:
        with gzip.open(paf_file, 'rb') as paf:
            collection = {}
            taxid_sum_dict = {}
            for _index, line in enumerate(paf):
                aln = Alignment(line)
                if aln.score >= score and int(aln.query_length) >= min_read_length:
                    if aln.query_name not in collection.keys():
                        collection[aln.query_name] = {"taxid": aln.taxid,
                                                    "read_length": aln.query_length,
                                                    "ref_accession": aln.ref_accession,
                                                    "read_coverage": aln.read_coverage,
                                                    "aligned_coverage": aln.aligned_coverage,
                                                    "identity": aln.aligned_identity,
                                                    "score": aln.score}
                    elif collection[aln.query_name]['score'] < aln.score:
                        collection[aln.query_name] = {"taxid": aln.taxid,
                                                    "read_length": aln.query_length,
                                                    "ref_accession": aln.ref_accession,
                                                    "read_coverage": aln.read_coverage,
                                                    "aligned_coverage": aln.aligned_coverage,
                                                    "identity": aln.aligned_identity,
                                                    "score": aln.score}
                # if _index == 200:
                #     break

            for k, _ in collection.items():
                if collection[k]['taxid'] not in taxid_sum_dict.keys():
                    taxid_sum_dict[collection[k]['taxid']] = {"count": 1}
                else:
                    taxid_sum_dict[collection[k]['taxid']]['count'] += 1
            
            taxid_sum_dict_sorted = {taxid: count for taxid, count in sorted(
                taxid_sum_dict.items(), reverse=True, key=lambda item: item[1]['count'])}

            _sum_taxonomy = []
            _metaphlan_like = []
            for k, v in taxid_sum_dict_sorted.items():
                _taxon_group = f"{k},{taxid_sum_dict_sorted[k]['count']}"
                get_taxonomy = ",".join(print_metaphlan_like_report(k).values())
                _sum_taxonomy.append(f"{_taxon_group},{get_taxonomy}")

                _taxon_group_metaphlan = f"{k}\t{taxid_sum_dict_sorted[k]['count']}"
                _get_taxonomy_metaphlan = "|".join(
                    print_metaphlan_like_report(k).values())
                # Remove | at the end
                _get_taxonomy_metaphlan[:_get_taxonomy_metaphlan.rfind("|")]
                _metaphlan_like.append(f"{_get_taxonomy_metaphlan}\t{_taxon_group_metaphlan}\t")

            logger.info("Writing reads")
            write_taxonomy_read(f"{output}_reads.csv", collection)
            
            logger.info("Writing report")
            write_taxonomy_report(f"{output}_report.csv", _sum_taxonomy)

            logger.info("Writing metaphlan-like-report")
            write_metaphlan_like_report(
                f"{output}_metaphlan_report.csv", f"{output}_report.csv")

    except IOError:
        logger.error("Could not open paf file")
# Testing with Pytest


def fake_single_cleaned_paf_line():
    paf_line = b'''9b4cfa9d-ef93-4c7f-8b66-364729c4b5a8\t3462\t2167\t3421\t+\tkraken:taxid|562|NZ_CP083266.1\t1865\t146\t1435\t1221\t1295\t60\tNM: i: 74 ms: i: 2155       AS: i: 2132
    nn: i: 0  tp: A: P  cm: i: 168        s1: i: 1027       s2: i: 899        de: f: 0.0424     rl: i: 50 cg: Z: 31M1I91M10D108M1D21M1D67M2D6M3D92M1I74M2D3M2D128M1D15M2D18M3D27M1I48M1I6M1D2M1D1M2D21M1D65M2D17M1D3M1D12M2D159M1D41M1D56M1I2M1I81M1D53M'''
    return paf_line


def test_parsing_paf_line():
    _line = fake_single_cleaned_paf_line()
    _alignment = Alignment(_line)
    assert _alignment.query_name == '9b4cfa9d-ef93-4c7f-8b66-364729c4b5a8'
    assert _alignment.query_length == '3462'
    assert _alignment.ref_name == '562|NZ_CP083266.1'
    assert _alignment.taxid == '562'
    assert _alignment.aligned_coverage == 35.27
    assert _alignment.score == 51.34


def test_swap_rank_dict():
    _dict = {1: 'a', 2: 'b', 3: 'c'}
    _swap_dict = {'a': 1, 'b': 2, 'c': 3}
    assert swap_rank_dict(_dict) == _swap_dict


def test_print_metaphlan_like_report():
    _taxonomy = print_metaphlan_like_report(525)
    assert _taxonomy['family'] == 'f__Acetobacteraceae'
    assert _taxonomy['kingdom'] == 'k__Bacteria'
    assert _taxonomy['species'] == 's__Acidocella facilis'

# pytest paf2taxid.py

if __name__ == '__main__':
    main()
