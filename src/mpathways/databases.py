#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""databases.py: Contains gene set databases for pathway analysis."""

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any, Union
from mbf.externals.util import download_file
from mbf.genomics.genes import Genes
from mbf.genomes import EnsemblGenome
from pypipegraph import FileGeneratingJob, Job
from mbf.genomics.annotator import Annotator
from .util import write_cls, write_gct
from pandas import DataFrame
from abc import abstractmethod
import pandas as pd
import pypipegraph2 as ppg2
import subprocess


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class GMTCollection:  # ExternalDataBase
    def __init__(self, name: str, genome: EnsemblGenome):
        self.name = name
        self.species = genome.species
        self.genome = genome
        self.cache_dir = Path("cache") / "gmt" / self.name / self.species
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._gmt = self.cache_dir / "input.gmt"
        self.dependencies: List[Job] = []

    def __str__(self):
        return f"{self.name} {self.species}"

    def __repr__(self):
        return f"{self.name} {self.species}"

    @property
    def filename(self):
        return self._gmt.absolute()

    @abstractmethod
    def write(self):
        pass

    @property
    def ensembl_gmt(self):
        return self.cache_dir / "ensembl.gmt"

    def write_ensembl(self):
        outfile = self.ensembl_gmt
        gmt_file = self._gmt
        genome = self.genome

        def __dump(outfile):
            with gmt_file.open("r") as inp:
                with outfile.open("w") as outp:
                    for line in inp.readlines():
                        splits = line[:-1].split("\t")
                        outp.write(f"{splits[0]}\t{splits[1]}")
                        for gene_name in splits[2:]:
                            stable_ids = list(genome.name_to_gene_ids(gene_name))
                            stable = "\t".join(stable_ids)
                            outp.write(f"\t{stable}")
                        outp.write("\n")

        return (
            ppg2.FileGeneratingJob(outfile, __dump)
            .depends_on(self.genome.download_genome())
            .depends_on(self.write())
            .depends_on(self.dependencies)  #
        )

    def get_depenbdencies(self) -> List[Job]:
        return self.dependencies


class GMTCollectionFromList(GMTCollection):  # ExternalDataBase
    def __init__(self, collections: List[GMTCollection], genome: EnsemblGenome, name: str = None):
        if name is None:
            name = "_".join([x.name for x in collections])
        super().__init__(name, genome)
        self.collections = collections
        for col in self.collections:
            self.dependencies.extend(col.dependencies)
            self.dependencies.append(col.write())

    def write(self):
        gmt = self._gmt
        filenames = [str(collection.filename) for collection in self.collections]

        def __dump(gmt):
            with gmt.open("w") as outp:
                for filename in filenames:
                    cmd = ["cat", f"{filename}"]
                    subprocess.check_call(cmd, stdout=outp)

        return ppg2.FileGeneratingJob(gmt, __dump).depends_on(self.dependencies)


class MSigDBCollection(GMTCollection):
    def __init__(
        self,
        collection_name: str,
        genome: EnsemblGenome,
        version: str = "7.1",
        subset: str = "all",
    ):
        name = ".".join([collection_name, subset, "v" + version])
        super().__init__(name, genome)
        if int(self.genome.revision) < 97:
            raise ValueError("Please use an Ensembl Genome from revision 97 onward.")
        try:
            v = float(version)
        except ValueError:
            raise ValueError(
                f"Cannot understand version {version}. It should be something like 7.1."
            )
        if v < 7:
            raise ValueError(
                "MSigDB Ensembl mapping only works for version 7.0 onward. Version was {version}."
            )
        self.version = version
        self.subset = subset
        self.collection_name = name
        self.input_file = self.cache_dir / (self.name + ".gmt")
        self.dependencies = [
            ppg2.ParameterInvariant(self.name, [self.collection_name, self.version, self.subset])
        ]

    def get_set_from_url(self, *_):
        url = f"https://data.broadinstitute.org/gsea-msigdb/msigdb/release/{self.version}/{self.name}.symbols.gmt"
        download_file(url, self._gmt.open("wb"))

    def write(self):
        return ppg2.FileGeneratingJob(self._gmt, self.get_set_from_url)


class IPACollection(GMTCollection):
    def __init__(
        self,
        name: str,
        genome,
    ):
        super().__init__(name, genome)
        if self.name == "ipa":
            self.url = "https://github.com/MarcoMernberger/ipa/raw/master/Functional%20Annotation%20IPA.txt"
        elif self.name == "ipa_reg":
            self.url = (
                "https://github.com/MarcoMernberger/ipa/raw/master/Regulator%20List%20IPA.txt"
            )
        else:
            raise ValueError(f"Cannot interpret name {name}.")
        self.cache_dir = Path("cache") / "gmt" / self.name / genome.species
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.input_file = self.cache_dir / f"{self.name}.txt"
        self.dependencies = [
            ppg2.ParameterInvariant(self.name, [self.name, str(self.input_file)]),
            self.load(),
        ]
        if self.name == "ipa":
            self.format_function = self.__format_ipa
        else:
            self.format_function = self.__format_ipa_reg

    def load(self) -> FileGeneratingJob:
        """
        Returns a FileGeneratingJob that downloads the Ensembl Chip annotation
        from MSigDB.

        The chip file is downloaded as input.chip in the cache directory.

        Returns
        -------
        FileGeneratingJob
            The job that creates the file.
        """

        url = self.url

        def __dump(input_file):
            download_file(url, input_file.open("wb"))

        return ppg2.FileGeneratingJob(self.input_file, __dump).depends_on(self.dependencies)

    @staticmethod
    def __format_ipa(line, name):
        splits = line.split("\t")
        genes = "\t".join(splits[1].split(", "))
        return "\t".join(
            [
                f"{splits[0]}",
                f"{name},{','.join([splits[2], splits[3], splits[4]])}",
                f"{genes}",
            ]
        )

    @staticmethod
    def __format_ipa_reg(line, name):
        splits = line.split("\t")
        genes = "\t".join(splits[2].split(", "))
        return "\t".join(
            [
                f"{splits[0]},{splits[1]}",
                f"{name},{','.join([splits[3], splits[4], splits[5]])}",
                f"{genes}",
            ]
        )

    def write(self) -> FileGeneratingJob:
        input_file = self.input_file
        format_fnc = self.format_function
        name = self.name.upper()

        def __dump(gmt):
            with gmt.open("w") as out:
                with input_file.open("r") as imp:
                    for line in imp.readlines()[1:]:
                        o = format_fnc(line[:-1], name)
                        out.write(o)
                        out.write("\n")

        return ppg2.FileGeneratingJob(self._gmt, __dump).depends_on(self.dependencies)


class MSigChipEnsembl:
    def __init__(self, species: str = "Homo_sapiens", version: str = "7.1"):
        """
        An ensembl chip object that takes care of file download and input.chip
        generation for GSEA.

        Parameters
        ----------
        version : str, optional
            The MSigDB version, by default "7.1"
        species : str, optional
            The species, by default "Homo_sapiens". Currently only supports
            Human, Mouse and Rat data.

        Raises
        ------
        ValueError
            If an unsupported species is provided.
        """
        if species == "Homo_sapiens":
            self.species = "Human"
            self.url = f"https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/{self.species}_ENSEMBL_Gene_MSigDB.v{version}.chip"
        elif species == "Mus_musculus":
            self.species = "Mouse"
            if version == "7.0":
                self.url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/Mouse_ENSEMBL_Gene_ID_MSigDB.v7.0.chip"
            elif version == "7.1":
                self.url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/Mouse_ENSEMBL_Gene_ID_to_Human_Orthologs_MSigDB.v7.1.chip"
            elif version == "7.2":
                self.url = f"https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/Mouse_ENSEMBL_Gene_ID_Human_Orthologs_MSigDB.v{version}.chip"
        elif species == "Rattus_norvegicus":
            self.species = "Rat"
            if version == "7.0":
                self.url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/Rat_ENSEMBL_Gene_ID_MSigDB.v7.0.chip"
            elif version == "7.1":
                self.url = "https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/Rat_ENSEMBL_Gene_ID_to_Human_Orthologs_MSigDB.v7.1.chip	30-Mar-2020 16:56"
            elif version == "7.2":
                self.url = f"https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/Rat_ENSEMBL_Gene_ID_Human_Orthologs_MSigDB.v{version}.chip"
        else:
            raise ValueError(
                f"Currently the species {species} is not supported. Check MsigDB chip files at https://data.broadinstitute.org/gsea-msigdb/msigdb/annotations_versioned/."
            )
        self.name = self.__class__.generate_name(species, version)
        self.version = version
        self.cache_dir = Path("cache") / "chip" / self.name
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._chip = self.cache_dir / "input.chip"
        self.dependencies = [ppg2.ParameterInvariant(self.name, [self.species, self.version])]

    @property
    def filename(self):
        return self._chip.absolute()

    @classmethod
    def generate_name(cls, species: str, version: str) -> str:
        return f"MSigChipEnsembl_{species}_{version}"

    def write(self) -> FileGeneratingJob:
        """
        Returns a FileGeneratingJob that downloads the Ensembl Chip annotation
        from MSigDB.

        The chip file is downloaded as input.chip in the cache directory.

        Returns
        -------
        FileGeneratingJob
            The job that creates the file.
        """
        outfile = self._chip
        url = self.url

        def __dump(outfile):
            download_file(url, outfile.open("wb"))

        return ppg2.FileGeneratingJob(outfile, __dump).depends_on(self.dependencies)


class CLSWriter:
    def __init__(
        self, name: str, phenotypes: Tuple[str, str], columns_a_b: Tuple[List[str], List[str]]
    ):
        self.name = name
        self.cache_dir = Path("cache") / "cls" / self.name
        self.columns_a_b = columns_a_b
        self.phenotypes = phenotypes
        self.dependencies = [
            ppg2.ParameterInvariant(self.name, list(phenotypes) + columns_a_b[0] + columns_a_b[1])
        ]
        self._cls = self.cache_dir / "input.cls"

    @property
    def filename(self):
        return self._cls.absolute()

    def write(self):
        return write_cls(self.cache_dir, self.phenotypes, self.columns_a_b, self.dependencies)


class GCTWriter:
    def __init__(
        self,
        comparison_name: str,
        genes_or_dataframe: Union[Genes, DataFrame],
        phenotypes: Tuple[str, str],
        columns_a_b: Tuple[List[str], List[str]],
        name: str = "Gct_df_default",
        dependencies: List[Job] = [],
    ):
        self.dependencies = dependencies
        if isinstance(genes_or_dataframe, Genes):
            self.name = f"Gct_{comparison_name}_{phenotypes[0]}_vs_{phenotypes[1]}"
            self.dependencies.append(genes_or_dataframe.load())
        else:
            self.name = f"{name}_{phenotypes[0]}_vs_{phenotypes[1]}"
        self.cache_dir = Path("cache") / "gct" / self.name
        self.columns_a_b = columns_a_b
        self.phenotypes = phenotypes
        self.genes_or_dataframe = genes_or_dataframe
        self.dependencies.append(
            ppg2.ParameterInvariant(
                self.name,
                list(phenotypes) + columns_a_b[0] + columns_a_b[1] + [self.name],
            )
        )
        self._gct = self.cache_dir / "input.gct"

    @property
    def filename(self):
        return self._gct.absolute()

    def write(self):
        return write_gct(
            self.genes_or_dataframe,
            self.cache_dir,
            self.phenotypes,
            self.columns_a_b,
            self.dependencies,
        )


def generate_collection(name: str, genome: EnsemblGenome) -> GMTCollection:
    if (
        name.startswith("h")
        or name.startswith("c1")
        or name.startswith("c2")
        or name.startswith("c3")
        or name.startswith("c4")
        or name.startswith("c5")
        or name.startswith("c6")
        or name.startswith("c7")
    ):
        if "." in name:
            splits = name.split(".")
            if name[-1].isdigit():
                return MSigDBCollection(
                    splits[0],
                    genome,
                    version=".".join([splits[-2][1:], splits[-1]]),
                    subset=".".join(splits[1:-2]),
                )
            else:
                return MSigDBCollection(splits[0], genome, subset=".".join(splits[1]))
        else:
            return MSigDBCollection(name, genome)
    elif name.startswith("ipa"):
        return IPACollection(name, genome)
    else:
        raise NotImplementedError(f"Don't know how to interpret this name: {name}.")


def interpret_collection(
    collection: Union[GMTCollection, List[GMTCollection], str, List[str]], genome: EnsemblGenome
) -> GMTCollection:
    """
    Turns a list of strings, GMTCollections or a single string into a single
    GMTCollection on which all subsequent methods can work.

    Parameters
    ----------
    collection : Union[GMTCollection, List[GMTCollection], str, List[str]]
        [description]
    genome : EnsemblGenome
        [description]

    Returns
    -------
    GMTCollection
        [description]

    Raises
    ------
    ValueError
        [description]
    """
    if isinstance(collection, str):
        collection = generate_collection(collection, genome)
    elif isinstance(collection, list):
        if isinstance(collection[0], GMTCollection):
            collection = GMTCollectionFromList(collection, genome)
        elif isinstance(collection[0], str):
            collection = GMTCollectionFromList(
                [generate_collection(x, genome) for x in collection],
                genome,
            )
        else:
            raise ValueError(f"Don't know how to interpret this collection: {collection}.")
    print("interpret collection", type(collection), collection)
    return collection
