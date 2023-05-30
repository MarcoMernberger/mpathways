#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""gsea.py: Contains wrapper for Gene Set Enrichment Analysis."""

from pathlib import Path
from typing import Optional, List, Dict, Tuple, Any, Union, Callable
from mbf.externals.util import download_zip_and_turn_into_tar_gzip, download_file
from mbf.genomics.genes import Genes
from mbf.genomics.annotator import Annotator
from mbf.genomes import EnsemblGenome
from .databases import (
    GMTCollection,
    MSigChipEnsembl,
    CLSWriter,
    GCTWriter,
    GMTCollectionFromList,
    interpret_collection,
)
from datetime import datetime
from pypipegraph2 import Job
from bs4 import BeautifulSoup
import pypipegraph2 as ppg2
import warnings
import os
import subprocess

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


global_instances: Dict[str, Any] = {}


class GSEA:
    def __init__(self):
        self.memory_in_mb = 16 * 1024
        self.collapse_values = ["No_Collapse", "Collapse", "Remap_Only"]
        self.collapse_modes = [
            "Max_probe",
            "Median_of_probes",
            "Sum_of_probes",
            "Mean_of_probes",
        ]
        self.weights = ["classic", "weighted", "weighted_p1.5", "weighted_p2"]
        self.metrics = [
            "Signal2Noise",
            "tTest",
            "Ratio_of_Classes",
            "Diff_of_Classes",
            "log2_Ratio_of_Classes",
            "Cosine",
            "Euclidean",
            "Manhattan",
            "Pearson",
        ]

    # @property
    # def java_path(self):
    #     path = list((self.path / f"GSEA_Linux_{self.version}").glob("jdk*"))[0] / "bin"
    #     return path

    # def unpacked_path(self):
    #     return self.store.get_unpacked_path("GSEA", self.version)

    # def unpack(self):
    #     self.store.unpack_version(self.name, self.version)
    #     self.__chmod_java()

    # def __chmod_java(self):
    #     for filename in self.java_path.iterdir():
    #         os.chmod(filename, 111)

    # def fetch_version(self, version: str, target_filename: Path) -> None:
    #     """
    #     Takes care of the tool download.

    #     Overrides the ExternalAlgorithm methood. Downloads the external method
    #     to the prebuild location specified by the corresponding
    #     ExternalAlgorithmStore and packs it into a tar.gz file.

    #     Parameters
    #     ----------
    #     version : str
    #         The tool version to be used.
    #     target_filename : Path
    #         The path to the local tar.gz file.
    #     """
    #     major = version[: version.rfind(".")]
    #     url = f"https://data.broadinstitute.org/gsea-msigdb/gsea/software/desktop/{major}/GSEA_Linux_{version}.zip"
    #     download_zip_and_turn_into_tar_gzip(
    #         url, target_filename, chmod_x_files=[f"GSEA_Linux_{version}/gsea-cli.sh"]
    #     )

    @property
    def name(self) -> str:
        """
        Returns the name of the encapsulated method for version handling.

        Returns
        -------
        str
            Name of the external method.
        """
        return "GSEA"

    def gsea_cmd(
        self,
        arguments: Union[None, List[str]],
    ):
        """
        Returns a command as a list of strings to be passed to subprocess.

        Constructs an executable command to be passed to subprocess from the
        output directory, the number of cores to be used and a list of command
        options. Overrides the ExternalAlgorithm method.

        Parameters
        ----------
        output_directory : Optional[Path]
            Path of output directory.
        ncores : int
            Number of cores to be used.
        arguments : List[str]
            List of string arguments for the command call.

        Returns
        -------
        List[str]
            Command to subprocess as a list of strings.
        """
        if arguments is None:
            arguments = []
        return [
            "gsea-cli.sh",
            "GSEA",
        ] + arguments

    def get_latest_version(self):
        """Returns the latest_version attribute."""
        return self.latest_version

    def verify_kwargs(self, **kwargs):
        accepted_parameters = [
            "norm",
            "permutations",
            "set_min",
            "set_max",
            "balance_rnd",
            "weight",
            "metric",
            "descending",
            "absolute_sorting",
            "top_x",
            "result_dir",
            "permutation_type",
            "median",
        ]
        warnings.simplefilter("default", UserWarning)
        for key in kwargs:
            if key not in accepted_parameters:
                warnings.warn(
                    f"Keyword {key} argument not recognized, it will be ignored. Accepted arguments are: {accepted_parameters}."
                )

    def check_arguments(
        self,
        collection: GMTCollection,
        gct: GCTWriter,
        chip: MSigChipEnsembl,
        clss: CLSWriter,
        result_dir: Path,
        rpt_label: str,
        **kwargs,
    ):
        arguments = []
        """
        collapse = kwargs.get("collapse", "No_Collapse")
        if collapse not in self.collapse_values:
            raise ValueError(f"Keyword argument 'collapse' mut be any of {self.collapse_values}, was {collapse}.")
        collapse_mode = kwargs.get("collapse_mode", "Max_probe")
        if collapse_mode not in self.collapse_modes:
            raise ValueError(f"Keyword argument 'collapse_mode' mut be any of {self.collapse_modes}, was {collapse_mode}.")
        """
        self.verify_kwargs(**kwargs)
        collapse = "Collapse"
        collapse_mode = "Max_probe"
        do_normalize = kwargs.get("norm", True)
        norm = "meandiv" if do_normalize else "None"
        permutations = kwargs.get("permutations", 1000)
        perm_type_auto = (
            "phenotype"
            if ((len(clss.columns_a_b[0]) >= 7) and (len(clss.columns_a_b) >= 7))
            else "gene_set"
        )
        perm_type = kwargs.get("permutation_type", perm_type_auto)
        set_min = kwargs.get("set_min", 5)
        set_max = kwargs.get("set_max", 3000)
        rnd = kwargs.get("balance_rnd", False)
        rnd_type = "equalize_and_balance" if rnd else "no_balance"
        weight = kwargs.get("weight", "weighted")
        if weight not in self.weights:
            raise ValueError(
                f"Keyword argument 'weight' mut be any of {self.weights}, was {weight}."
            )
        metric = kwargs.get("metric", "Signal2Noise")
        if metric not in self.metrics:
            raise ValueError(
                f"Keyword argument 'metric' mut be any of {self.metrics}, was {metric}."
            )
        descending = "descending" if kwargs.get("descending", True) else "ascending"
        sorting = "real" if not kwargs.get("absolute_sorting", False) else "abs"
        topx = kwargs.get("top_x", 200)
        median = kwargs.get("median", False)
        arguments = [
            "-res",
            str(gct.filename),  # this needs to come from genes
            "-cls",
            str(clss.filename),  # this needs to come from comparison
            "-chip",
            str(chip.filename),  # we try using the MSigDB chips chip.chip
            "-gmx",
            str(collection.filename),  # this needs to be done once for the GMTset collection.gmt
            "-set_max",
            str(set_max),
            "-set_min",
            str(set_min),
            "-out",
            str(result_dir),
            "-collapse",
            collapse,
            "-mode",
            collapse_mode,
            "-norm",
            norm,
            "-nperm",
            str(permutations),
            "-permute",
            perm_type,
            "-rnd_type",
            rnd_type,
            "-scoring_scheme",
            weight,
            "-rpt_label",
            rpt_label,
            "-metric",
            metric,
            "-sort",
            sorting,
            "-order",
            descending,
            "-create_svgs",
            "true",
            "-create_gcts",
            "true",
            "-include_only_symbols",
            "true",
            "-make_sets",
            "true",
            "-median",
            str(median).lower(),
            "-num",
            "100",
            "-plot_top_x",
            str(topx),
            "-rnd_seed",
            "timestamp",
            "-save_rnd_lists",
            "false",
            "-zip_report",
            "true",
        ]
        return arguments

    def get_file_writer(
        self,
        comparison_name: str,
        genes: Genes,
        phenotypes: Tuple[str, str],
        columns_a_b: Tuple[List[str], List[str]],
        chip: Union[MSigChipEnsembl, None],
        dependencies: List[Job],
    ):
        cls_writer = CLSWriter(comparison_name, phenotypes, columns_a_b)
        gct_name = f"Gct_{genes.name}_{phenotypes[0]}_vs_{phenotypes[1]}"
        gct_writer = GCTWriter(
            comparison_name, genes, phenotypes, columns_a_b, gct_name, dependencies
        )
        return cls_writer, gct_writer

    def __clean_date_folder(self, outdir, now, rpt_label):
        def __clean():
            month = now.strftime("%b").lower()
            for subdir in outdir.iterdir():
                if subdir.name.startswith(month) and subdir.is_dir():
                    subdir.rmdir()
                if str(subdir.name).startswith(rpt_label):
                    new_name = str(subdir.parent / rpt_label)
                    subdir = subdir.rename(new_name)
                    self.__create_main_index(Path(new_name))

        return __clean

    def __get_index_file(self, outdir, rpt_label):
        for subdir in outdir.iterdir():
            if subdir.startswith(rpt_label):
                return subdir / "index.html"
        return None

    def __create_main_index(self, folder: Path):
        index_file = folder / "index.html"
        with (folder.parent / "index.html").open("w") as outp:
            soup = BeautifulSoup(index_file.open("r"), "html")
            for x in soup.find_all("a"):
                href = x.attrs["href"]
                if not href.startswith("http"):
                    x.attrs["href"] = str(Path(folder.name) / href)
            outp.write(soup.prettify())

    def run_on_counts(
        self,
        genes: Genes,
        comparison_name: str,
        phenotypes: Tuple[str, str],
        columns_a_b: Tuple[List[str], List[str]],
        collection: Union[GMTCollection, List[GMTCollection], str, List[str]],
        genome: EnsemblGenome,
        chip: Optional[MSigChipEnsembl] = None,
        annotators: List[Annotator] = [],
        **kwargs,
    ):
        dependencies = []
        if chip is None:
            chip = MSigChipEnsembl(genome.species, "7.1")
        for anno in annotators:
            dependencies.append(genes.add_annotator(anno))
        clss, gct = self.get_file_writer(
            comparison_name, genes, phenotypes, columns_a_b, chip, dependencies
        )
        collection = interpret_collection(collection, genome)
        result_dir = kwargs.get(
            "result_dir",
            Path("results/GSEA") / genes.name / comparison_name / collection.name,
        )
        result_dir = result_dir.absolute()
        result_dir.mkdir(parents=True, exist_ok=True)
        now = datetime.now()
        rpt_label = now.strftime("%Y%m%d%H%M")
        index_html = (result_dir / "index.html").relative_to("/project")
        arguments = self.check_arguments(
            collection, gct, chip, clss, result_dir, rpt_label, **kwargs
        )
        job = self.run_job(
            arguments,
            result_dir,
            call_afterwards=self.__clean_date_folder(result_dir, now, rpt_label),
            # additional_files_created=[index_html],
        )
        job.depends_on(chip.write())
        job.depends_on(collection.write())
        job.depends_on(clss.write())
        job.depends_on(gct.write())
        return job, index_html

    def run_job(self, arguments: List[str], result_dir: Path, call_afterwards: Callable):
        outfile = result_dir / "index.html"
        cmd = self.gsea_cmd(arguments)
        sentinel = result_dir / "sentinel.txt"
        sentinel = sentinel.relative_to("/project")

        def __call(sentinel: Path):
            try:
                with sentinel.open("w") as outp:
                    subprocess.run(cmd, capture_output=True, text=True, check=True)
                    call_afterwards()
                    outp.write("GSEA command:\n" + " ".join(cmd))

            except subprocess.CalledProcessError:
                print(" ".join(cmd))
                raise

        return ppg2.FileGeneratingJob(sentinel, __call)

    def __repr__(self) -> str:
        return f"GSEA"
