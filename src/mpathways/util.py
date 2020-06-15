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