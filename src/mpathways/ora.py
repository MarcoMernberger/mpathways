#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""ora.py: Contains tools for Over-Representation Analysis."""

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union
from mbf_genomics.genes import Genes
from mbf_genomes import EnsemblGenome
from .databases import GMTCollection, interpret_collection
from .util import fdr_control_benjamini_hochberg
from datetime import datetime
from pypipegraph import Job, FileGeneratingJob, CachedAttributeLoadingJob
from mplots import MPPlotJob
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pypipegraph as ppg
import scipy
import numpy as np
import math

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class ORAHyper:
    """
    Class that offers an ORA analysis using hypergeometrix enrichment.

    The class takes care of gene set preparation and prerequisites exactly
    once and allows for analysis of multiple DE gene sets.
    Therefore, an instance of HE is defined by the genome, background gene
    set and the collection used for analysis.

    Parameters
    ----------
    name : str
        Name of the HE instance.
    genome : EnsemblGenome
        The Ensembl genome used.
    background_genes : Genes
        Genes instance with all measured genes, to be used as background.
    collection : Union[GMTCollection, List[GMTCollection], str, List[str]]
        Gene set database to be used. If a list of GMTCollections is given,
        it is combined to a single GMTCollection. If a str or list of str
        is given, str are interpreted as names and resolved if possible.
    """

    def __init__(
        self,
        name: str,
        genome: EnsemblGenome,
        background_genes: Genes,
        collection: Union[GMTCollection, List[GMTCollection], str, List[str]],
    ) -> None:
        self.name = name
        self.genome = genome
        self.collection = interpret_collection(collection, genome)
        self.background_genes = background_genes
        self.all_sets = set()
        self.dependencies = [
            self.collection.write_ensembl(),
            self.background_genes.load(),
        ]
        self.cache_dir = Path("cache") / "HE" / self.name
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def trim_sets(self) -> CachedAttributeLoadingJob:
        """
        Trims the collections and the background gene sets to the genes that
        could have benn observed as they were actually measured.

        A triple of 3 objets is loaded as attribute called sets_allinset_desc,
        consisting of sets (a dictionary of trimmed sets to use), all_in_set
        (a set of all observed genes) and desc (a dictionary of set descriptions).

        Returns
        -------
        CachedAttributeLoadingJob
            Job that loads the attribte.
        """

        def calc():
            background = set(self.background_genes.df["gene_stable_id"])
            sets = {}
            allinset = set()
            desc = {}
            with self.collection.ensembl_gmt.open("r") as inp:
                for line in inp.readlines():
                    splits = line[:-1].split("\t")
                    ids = set(splits[2:]).intersection(background)
                    allinset.update(ids)
                    sets[splits[0]] = ids
                    desc[splits[0]] = splits[1]
            return (sets, allinset, desc)

        return ppg.CachedAttributeLoadingJob(
            self.cache_dir / "sets", self, "sets_allinset_desc", calc
        ).depends_on(self.dependencies)

    def run(self, genes: Genes, **kwargs) -> FileGeneratingJob:
        """
        Runs the analysis and creates an output table for the results.

        Returns a job that creates the table.

        Parameters
        ----------
        genes : Genes
            The selected genes to use.

        Returns
        -------
        FileGeneratingJob
            Jov that generates the file.
        """
        outdir = kwargs.get("result_dir", (genes.result_dir / "ORA"))
        outdir.mkdir(parents=True, exist_ok=True)
        outfile = outdir / f"{genes.name}_{self.name}.tsv"

        def __run():
            def get_gene_name(x):
                try:
                    return genes.genome.genes[x].name
                except IndexError:
                    return "<no name>"

            result = {
                "Group": [],  # whan functional gene groups collection are we talking about?
                "Set": [],  # and what was the name of the set
                "Set size": [],  # how many are in this set
                "Background size": [],  # and how many are there in total (in set + not in set)
                "Observed overlap": [],  # overlap query ^ this set
                "Query size": [],  # overlap query ^ background
                "Query size including non-annotated": [],
                "p-value": [],
                "Description": [],  # where do you get more information about this set?
                "Genes": [],
                "GeneNames": [],
            }
            sets_to_test, possiblyseen, description = self.sets_allinset_desc
            df = genes.df
            query_set = set(df["gene_stable_id"])
            query_size_including_non_annotated = len(query_set)
            print("len query set before", len(query_set))
            query_set = query_set.intersection(
                possiblyseen
            )  # because only these could have possibly hit
            no_of_trials = len(query_set)
            print("len query set after", no_of_trials)
            total_gene_count = len(possiblyseen)
            for set_name in sets_to_test:
                geneset = sets_to_test[set_name]  # background is already intersected
                drawn = geneset.intersection(query_set)
                no_of_white_balls_in_urn = len(geneset)
                no_of_white_balls_drawn = len(drawn)
                no_of_black_balls_in_urn = total_gene_count - no_of_white_balls_in_urn
                p_value = scipy.stats.hypergeom(
                    M=no_of_white_balls_in_urn
                    + no_of_black_balls_in_urn,  # total number of balls
                    n=no_of_white_balls_in_urn,  # number of white balls
                    N=no_of_trials,  # no of balls drawn
                ).sf(no_of_white_balls_drawn - 1)
                result["Group"].append(self.collection.name)
                result["Set"].append(set_name)
                result["Set size"].append(no_of_white_balls_in_urn)
                result["Background size"].append(total_gene_count)
                result["Observed overlap"].append(no_of_white_balls_drawn)
                result["Query size"].append(no_of_trials)
                result["Query size including non-annotated"].append(
                    query_size_including_non_annotated
                )
                result["p-value"].append(p_value)
                result["Description"].append(description[set_name])
                genes_str = ", ".join(list(sorted(drawn)))
                result["Genes"].append(genes_str)
                gene_names = [get_gene_name(x) for x in list(sorted(drawn))]
                gene_names = ", ".join(gene_names)
                result["GeneNames"].append(gene_names)
            df = pd.DataFrame(result)
            fdr_control_benjamini_hochberg(df, "p-value", "Benjamini")
            partial_order = ["Set", "Group", "Benjamini", "Description"]
            column_order = partial_order + [
                x for x in df.columns if x not in partial_order
            ]
            df = df[column_order]
            df = df.sort_values("Benjamini", ascending=True)
            df = df[df["Observed overlap"] > 0]
            df.to_csv(str(outfile), sep="\t", index=False)

        return (
            ppg.FileGeneratingJob(outfile, __run)
            .depends_on(self.trim_sets())
            .depends_on(self.dependencies)
            .depends_on(genes.load())
        )

    @classmethod
    def plot_bars(
        cls, ora_job: FileGeneratingJob, topx: int = 50, dependencies=[], **kwargs
    ) -> MPPlotJob:
        """
        Plot a bar plot of enrichment terms on the x-axis and log-p on the y axis.

        Parameters
        ----------
        ora_job : FileGeneratingJob
            Job that creates the output table which is read for the plot generation.
        topx : int, optional
            Plot the top 50 significant hits, by default 50.
        dependencies : list, optional
            Additional dependencies, by default [].

        Returns
        -------
        MPPlotJob
            The Job that creates the file.
        """
        ora_result_file = ora_job.job_id
        plot_filename = Path(ora_result_file).with_suffix(f".top{topx}.png")
        dependencies.append(ora_job)

        def calc():
            df = pd.read_csv(ora_result_file, sep="\t")
            df = df[df["Benjamini"] <= 0.05].sort_values("Benjamini", ascending=True)
            df = df[: min(len(df), topx)]
            return df

        def plot(df):
            figsize = kwargs.get("figsize", (10, 10))
            title = kwargs.get("title", plot_filename.stem)
            fig = plt.figure(figsize=figsize)
            if len(df):
                df["corr p_value"] = -1 * np.log10(df["Benjamini"].values)
                df = df.sort_values("corr p_value", ascending=True)
                df["Set"] = pd.Categorical(
                    df["Set"].values, list(set(df["Set"].values))
                )
                y = range(len(df))
                plt.barh(
                    y,
                    df["corr p_value"].values,
                    height=0.8,
                    tick_label=df["Set"].values,
                )
                lims = plt.ylim()
                plt.gca().axvline(
                    [-1 * math.log(0.05, 10)], ymin=lims[0], ymax=lims[1], color="k"
                )
                plt.xlabel(r"$-log_{10}$ (corrected p-value)")
                plt.ylabel("enriched gene sets")
                plt.title(title)
            else:
                plt.title("no significant entries")
            plt.tight_layout()
            return fig

        return MPPlotJob(
            plot_filename,
            calc,
            plot,
            calc_args=[topx],
            plot_args=sorted(list(kwargs.items())),
        ).depends_on(dependencies)

