#!/usr/bin/env python3
import argparse
import os
import sys
import scipy.stats as stats
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests


class Term(object):
    """docstring for Term."""

    def __init__(self, term, level=0):
        super(Term, self).__init__()
        self.term = term.upper()
        self.gene_associated = []
        # self.genes_in_set = 0           # number of genes in the analysis set for term
        # genes (id) belonging to this term in gene set
        self.genes_in_set = []
        # genes (id) associated with this term in gene universe
        self.genes_in_universe = []
        self.level = level              # hierarchy level

    def __str__(self):
        return "TERM: {}\nN ASSOCIATED GENES: {} | N GENES IN SET:{} | N GENES IN UNIVERSE: {} \nGENES IN UNIVERSE: \n{}".format(self.term, len(self.gene_associated), len(self.genes_in_set), len(self.genes_in_universe), self.genes_in_universe)

    def add_associated_gene(self, gene_id):
        if gene_id not in self.gene_associated:
            self.gene_associated.append(gene_id)

    def add_genes_in_set(self, gene_list):
        # count the number of genes belonging to this
        self.genes_in_set = list(set(self.gene_associated).intersection(
            gene_list))  # <- sort before would speed up?
        # self.genes_in_set = len(self.genes_in_set)

    def add_genes_in_universe(self, gene_universe_list):
        self.genes_in_universe = list(
            set(self.gene_associated).intersection(gene_universe_list))

    def enrich(self, total_genes_set, total_genes_univ, test="hypergeometric"):
        # total number of genes on genome/analyzed
        n_univ = total_genes_univ
        # total number of genes with term associated in universe/population
        n_assoc_univ = len(self.genes_in_universe)
        # total number of genes in the analyzed subset
        n_set = total_genes_set
        # total number of genes with associated term in the analyzed subset
        n_assoc_set = len(self.genes_in_set)

        if(test == "fisher"):
            contingency_table = [[n_assoc_set,  n_set - n_assoc_set],
                                 [n_assoc_univ - n_assoc_set,
                                  n_univ - n_set - n_assoc_univ - n_assoc_set]
                                 ]
            oddsratio, pvalue = stats.fisher_exact(
                contingency_table, alternative='less')

            return [self.term.lower(),
                    n_assoc_set, n_set, n_assoc_univ, n_univ, pvalue]

        if(test == "hypergeometric"):
            sf = stats.hypergeom.sf(
                n_assoc_set - 1, n_univ, n_assoc_univ, n_set)
            # return data
            return [self.term.lower(),
                    n_assoc_set, n_set, n_assoc_univ, n_univ, sf]


class TermList(object):
    """docstring for TermList."""

    def __init__(self, file_name):
        super(TermList, self).__init__()
        self.terms_list = {}
        self.read_term_association_file(file_name)
        self.result_table = pd.DataFrame(
            columns=["term", "genes_in_set", "total_genes_in_set",
                     "genes_in_universe",
                     "total_genes_in_universe", "pvalue"])

    def read_term_association_file(self, file_name):
        with open(file_name) as mapping_file:
            for line in mapping_file:
                if line.startswith('#'):
                    continue

                # split line
                cols = line.rstrip('\n').split('\t')
                # col[0]: gene_id col[1] term
                # split term and basic hierarchy
                term_cols = cols[1].split(":")
                term_k = term_cols[1].upper()

                if term_k not in self.terms_list:
                    # create the term with name (upper) & hierarchy level
                    self.terms_list[term_k] = Term(str(term_k), term_cols[0])
                    self.terms_list[term_k].add_associated_gene(
                        cols[0])  # add the associated gene
                else:
                    self.terms_list[term_k].add_associated_gene(
                        cols[0])  # add the associated gene to the term

    def __process_gene_set_and_universe(self, gene_set_list,
                                        gene_universe_list):
        for term_k in self.terms_list:
            self.terms_list[term_k].add_genes_in_set(gene_set_list)
            self.terms_list[term_k].add_genes_in_universe(gene_universe_list)

    def __enrich_terms(self, n_gene_set, n_gene_univ, test):
        for term_k in self.terms_list:
            term = self.terms_list[term_k]
            # this should return
            res = term.enrich(n_gene_set, n_gene_univ, test)
            self.result_table.loc[len(self.result_table)] = res

    def __adj_pvalues(self, method):
        # correct p-values
        mt_re = multipletests(
            self.result_table['pvalue'], method=method)
        # add
        self.result_table['adj_pvalue'] = pd.Series(
            mt_re[1], index=self.result_table.index)

    def enrichment_test(self, gene_set_list, gene_universe_list,
                        test, correction):
        self.__process_gene_set_and_universe(gene_set_list, gene_universe_list)
        print("Total number of terms: {}".format(
            len(self.terms_list)), file=sys.stderr)

        self.__enrich_terms(len(gene_set_list), len(gene_universe_list), test)
        self.__adj_pvalues(correction)

    def save_enrichment_result(self, filename):
        self.result_table.to_csv(filename, sep="\t",
                                 index=False)


def filter_value(filter_type, filter_value, value):
    if filter_type == "l":
        if value < filter_value:
            return False
        else:
            return True
    elif filter_type == "h":
        if value > filter_value:
            return False
        else:
            return True


arg_parser = argparse.ArgumentParser(
    description="Metabolic class enrichment on provided metabolic mapping\n(AKA a pretty generic statistical term enrichment tool with basic hierarchy support)")
arg_parser.add_argument("--gene_list", "-g",
                        help="File containing gene list (subset of gene universe)",
                        required=True)

arg_parser.add_argument("--gene_universe", "-u",
                        help="File containing all genes considered in the expression experiment",
                        required=True)

arg_parser.add_argument("--gene_terms_mapping", "-m",
                        help="File containing gene_id (column 1 in gene_list) and LEVEL:Term. Tab separated",
                        required=True)

arg_parser.add_argument("--statistical_test", "-t",
                        help="Statistical test used for determine significance.",
                        default="hypergeometric",
                        choices=["hypergeometric", "fisher"],
                        required=False)

arg_parser.add_argument("--multiple_comparison", "-c",
                        help="Multiple testing correction",
                        default="bonferroni",
                        choices=["bonferroni", "sidak", "holm-sidak", "holm",
                                 "simes-hochberg", "hommel", "fdr_bh",
                                 "fdr_by", "fdr_tsbh", "fdr_tsbky"],
                        required=False)

arg_parser.add_argument("--alpha_filtering", "-a",
                        help="""Specify a filtering option using one of this combination. Gene list is expected to have a second numeric column\n
                        t:n - Filter top n (percent of total) gene based on numeric value of column two from 'gene_list' file\n
                        c:n - Filter top n (number) gene based on numeric value of column two from 'gene_list' file
                        h:n - Filter all genes with values contained in column two higher than n\n
                        l:n - Filter all genes with values contained in column two lower than n (exclusion)""",
                        required=False)
arg_parser.add_argument("--print_associated_ids", "-i",
                        help="In addition to all the normal information, output a column containing all the associated gene id",
                        required=False,
                        action='store_true')
arg_parser.add_argument("--output_file", "-O",
                        help="Save result to file",
                        default=sys.stdout)

cmd_args = arg_parser.parse_args()

# read mappings
terms_mapping_obj = TermList(cmd_args.gene_terms_mapping)

filter_t = "NO"
filter_a = 0
if cmd_args.alpha_filtering is not None:
    filter_t = cmd_args.alpha_filtering.split(":")[0]
    filter_a = int(cmd_args.alpha_filtering.split(":")[1])

# read gene set
gene_set = []
total_gene_in_set_file = 0
with open(cmd_args.gene_list) as genes_file:
    # store everything in a "table"

    for line in genes_file:
        if line.startswith('#'):
            continue

        cols = line.rstrip("\n").split("\t")
        if not filter_t == "NO":
            if filter_t == "h" or filter_t == "l":
                # RIVEDERE!
                if not filter_value(filter_t, filter_a, float(cols[1])):
                    gene_set.append(cols[0])
                    total_gene_in_set_file += 1
            else:
                gene_set.append({'gene_id': cols[0], 'value': float(cols[1])})
                total_gene_in_set_file += 1
        else:
            # list is already ok
            gene_set.append(cols[0])
            total_gene_in_set_file += 1
    # if we are handling top genes
    if filter_t == "t" or filter_t == "c":
        print("Total genes in set: {}".format(
            total_gene_in_set_file), file=sys.stderr)
        # calculate the number of gene to get
        num_of_genes = 0
        if filter_t == "t":
            # percentage
            num_of_genes = int(
                round((total_gene_in_set_file * filter_a) / 100))
            print("Filter: Top {}% ({})".format(
                filter_a, num_of_genes), file=sys.stderr)
        elif filter_t == "c":
            # Top N
            num_of_genes = filter_a
            print("Filter: Top {}".format(filter_a), file=sys.stderr)
        # sort by value
        sorted_list = sorted(gene_set, key=lambda k: k['value'], reverse=True)
        # subset the gene set and get only keys
        gene_set = map(lambda x: x['gene_id'], sorted_list[0:num_of_genes])
        # change total number of genes
        total_gene_in_set_file = num_of_genes


print("Total genes in set: {}".format(total_gene_in_set_file), file=sys.stderr)
# read gene universe
gene_univ = []
total_gene_univ_file = 0
with open(cmd_args.gene_universe) as univ_file:
    for line in univ_file:
        if line.startswith('#'):
            continue

        cols = line.rstrip("\n").split("\t")
        gene_univ.append(cols[0])
        total_gene_univ_file += 1

print("Test: {}".format(cmd_args.statistical_test))
print("Multiple testing correction: {}".format(cmd_args.multiple_comparison))
print("Total genes in universe: {}".format(
    total_gene_univ_file), file=sys.stderr)
# set interesting genes and gene universe
terms_mapping_obj.enrichment_test(gene_set, gene_univ,
                                  cmd_args.statistical_test,
                                  cmd_args.multiple_comparison)
terms_mapping_obj.save_enrichment_result(cmd_args.output_file)
