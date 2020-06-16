#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""util.py: Contains some utility functions."""

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union
from pypipegraph import Job, FileGeneratingJob
from mbf_genomics.genes import Genes
from pandas import DataFrame
import pandas as pd
import pypipegraph as ppg
import mbf_genomics
import scipy
import numpy as np


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def write_cls(output_directory: Path, phenotypes: Tuple[str, str], columns_a_b: Tuple[List[str], List[str]], dependencies: List[Job] = []) -> FileGeneratingJob:
    """
    Creates a Job that writes class information file for GSEA at a specified
    folder. A file named imput.cls is created.

    Parameters
    ----------
    output_directory : Path
        The output directory in which an input.cls file is created.
    phenotypes : Tuple[str, str]
        The phenotype/class names of the groups to be compared.
    columns_a_b : Tuple[List[str], List[str]]
        The DataFrame columns of the relevant expression values.
    dependencies : List[Job], optional
        List of prerequisite jobs on which the , by default []

    Returns
    -------
    FileGeneratingJob
        The job that creates the file.
    """
    output_directory.mkdir(parents=True, exist_ok=True)
    outfile = output_directory / "input.cls"

    def __dump():
        with outfile.open("w") as handle:
            handle.write(f"{(len(columns_a_b[0]) + len(columns_a_b[1]))} 2 1\n")
            handle.write(f"# {phenotypes[0]} {phenotypes[1]}\n")
            handle.write(" ".join(["0"] * len(columns_a_b[0])))
            handle.write(" ")
            handle.write(" ".join(["1"] * len(columns_a_b[1])))
            handle.write("\n")
            handle.close()

    return ppg.FileGeneratingJob(outfile, __dump).depends_on(dependencies)


def write_gct(genes_or_dataframe: Union[Genes, DataFrame], output_directory: Path, phenotypes: Tuple[str, str], columns_a_b: Tuple[List[str], List[str]], dependencies: List[Job] = []) -> FileGeneratingJob:
    """
    Creates a Job that writes expression data for GSEA at a specified
    folder. A file named imput.gct is created.

    Parameters
    ----------
    output_directory : Path
        The output directory in which an input.gct file is created.
    phenotypes : Tuple[str, str]
        The phenotype/class names of the groups to be compared.
    columns_a_b : Tuple[List[str], List[str]]
        The DataFrame columns of the relevant expression values.
    dependencies : List[Job], optional
        List of prerequisite jobs on which the , by default []

    Returns
    -------
    FileGeneratingJob
        The job that creates the file.
    """
    output_directory.mkdir(parents=True, exist_ok=True)
    outfile = output_directory / "input.gct"
    if isinstance(genes_or_dataframe, Genes):
        dependencies.append(genes_or_dataframe.add_annotator(mbf_genomics.genes.annotators.Description()))
    def __dump():
        df = genes_or_dataframe
        if isinstance(genes_or_dataframe, Genes):
            df = genes_or_dataframe.df.copy()
        elif isinstance(genes_or_dataframe, DataFrame):
            df = df.copy()
        else:
            raise ValueError(f"Parameter genes_or_dataframe must be an instance of Genes or DataFrame, was {type(genes_or_dataframe)}.")
        with outfile.open("w") as handle:
            handle.write("#1.2\n")
            handle.write(f"{len(df)}\t{len(columns_a_b[0]) + len(columns_a_b[1])}\n")
            handle.write("ProbeName\tDescription\t")
            handle.write("\t".join(columns_a_b[0] + columns_a_b[1]))
            handle.write("\n")
            df = df.rename(columns={"gene_stable_id": "NAME"})
            description = [f"{x} {y}" for x, y in zip(df["name"], df["description"])]
            df["Description"] = description
            df = df[["NAME", "Description"] + columns_a_b[0] + columns_a_b[1]]
            df = df.fillna(0)
            for _, row in df.iterrows():
                handle.write("\t".join([str(x) for x in row]) + "\n")

    return ppg.FileGeneratingJob(outfile, __dump).depends_on(dependencies)


def _benjamini_hochberg(col):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    A clever implementation that manipulates the results so that
    anything below your significance threshold is significant - even if
    the original BH calculation had a part where the FDR rose.
    (Remember, we reject all null-hypotheses below the one with FDR < alpha
    with the *highest* p-value - even if their FDR is >= alpha!)

    This is a direct translation from the R code, in essence.
    """
    p = np.asfarray(col)
    if (pd.isnull(p).any()):
        orig_idx = np.array(list(range(len(p))))
        is_nan = pd.isnull(p)
        q_ommiting_nans = _benjamini_hochberg(p[~is_nan])
        indices_without_nan = orig_idx[~is_nan]
        result = np.empty(len(p))
        result[:] = np.nan
        for q, idx in zip(q_ommiting_nans, indices_without_nan):
            result[idx] = q
        return result
    else:
        by_descend = p.argsort()[::-1]
        by_orig = by_descend.argsort()
        steps = float(len(p)) / np.arange(len(p), 0, -1)
        q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
        return q[by_orig]

def fdr_control_benjamini_hochberg(df,
                                   p_value_column='p-value',
                                   output_column='benjamini_hochberg',
                                   drop_output_if_exists=False):
    """Calculate the benjamini hochberg adjusted p-values in order to control false discovery rate
    Adds two columns, output_column and output_column + '_rank', since you need to order
    the p-values later to decide if a given value is significant.
    """
    if output_column in df.columns:
        if drop_output_if_exists:
            df = df.drop(output_column, axis=1)
        else:
            raise ValueError(
                "Dataframe already has a column called %s" % output_column)
    col = df[p_value_column]
    bh = _benjamini_hochberg(col)
    df.loc[:, output_column] = bh

def hypergeometric_test(query_set, reference_set, background_set):
    """Query set is what you observed, reference set is what you compare against,
    background set is what you could have observed (genes on array...) - and which were annotated(!)"""
    query_set = query_set.intersection(background_set)
    reference_set = reference_set.intersection(background_set)

    drawn = reference_set.intersection(query_set)
    total_gene_count = len(background_set)

    no_of_trials = len(query_set)
    no_of_white_balls_in_urn = len(reference_set)
    no_of_white_balls_drawn = len(drawn)
    no_of_black_balls_in_urn = total_gene_count - no_of_white_balls_in_urn
    return scipy.stats.hypergeom(
            M=no_of_white_balls_in_urn + no_of_black_balls_in_urn, # total number of balls
            n = no_of_white_balls_in_urn, #number of white balls
            N = no_of_trials # no of balls drawn
            ).sf(
                    no_of_white_balls_drawn  #no of white balls drawn
                    -1)


def multi_hypergeom_test(genome, query_set, function_gene_groups_or_list_of_such = None, background_set = None):
    """Test a query set against multiple sets from functional.databases.
    Returns a pandas.DataFrame(group, set, benjamini, p, overlap_count,
    sorted by benjamini

    """
    if function_gene_groups_or_list_of_such is None:
            function_gene_groups_or_list_of_such = databases.get_default_groups()
    query_set = set(query_set)
    list_of_gene_groups = check_function_gene_groups_or_list_of_such(function_gene_groups_or_list_of_such)
    sets_by_func_group = {}
    all_annotated_genes = set()
    for func_group in list_of_gene_groups:
        sets_by_func_group[func_group.name] = func_group.get_sets(genome)
        for v in sets_by_func_group[func_group.name].values():
            all_annotated_genes.update(v)
    if background_set is not None:
        background_set = background_set.intersection(all_annotated_genes)
    else:
        background_set = all_annotated_genes
    query_set = query_set.intersection(background_set)
    result = {'group': [], 'set': [], 'p': [], 'overlap count': [], 'intersection': []}
    for func_group_name in sets_by_func_group:
        for set_name, set_genes in sets_by_func_group[func_group_name].items():
            set_genes = set(set_genes)
            result['group'].append(func_group_name)
            result['set'].append(set_name)
            p = hypergeometric_test(query_set, set_genes, background_set)
            if p > 1:
                raise ValueError("p > 1,. was %.15f" % p)

            result['p'].append(p)
            if result['p'][-1] == 0:
                raise ValueError()
            intersection = query_set.intersection(set_genes)
            result['overlap count'].append(len(intersection))
            result['intersection'].append(", ".join(list(sorted([get_gene_name(genome, x) for x in intersection]))))
    res =  pd.DataFrame(result, )
    statistics.fdr_control_benjamini_hochberg(res, 'p', 'benjamini')
    res = res[['group','set','benjamini', 'p','overlap count', 'intersection']].sort_values('benjamini')
    return res