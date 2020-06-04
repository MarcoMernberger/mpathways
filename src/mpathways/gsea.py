#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""gsea.py: Contains wrapper for Gene Set Enrichment Analysis."""

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any
from mbf_externals import ExternalAlgorithm, ExternalAlgorithmStore
from mbf_externals.util import download_zip_and_turn_into_tar_gzip, download_file
import pandas as pd
import pypipegraph as ppg


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


class GSEA(ExternalAlgorithm):

    def __init__(
        self,
        version: str = "_last_used",
        store: Optional[ExternalAlgorithmStore] = None,
        **kwargs
    ):
        super().__init__(version=version, store=store, **kwargs)
        self.memory_in_mb = 16 * 1024
        """
        self.n_permutation = n_permutation
        self.name = name
        self.gene_groups = gene_groups
        self.genome = genome
        ppg.assert_uniqueness_of_object(self)
        self.output_path = Path("results") / "GSEA" / self.name
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.cache_path = Path("cache") / "GSEA" / self.name
        self.cache_path.mkdir(parents=True, exist_ok=True)
        self.gct_file = self.output_path / "input.gct"
        self.cls_file = self.output_path / "input.cls"
        self.chp_file = self.output_path / "input.chip"
        self.gmt_file = self.output_path / "input.gmt"
        self.do_collapse = collapse
        self.plot_top_x = plot_top_x
        self.set_min = set_min
        self.set_max = set_max
        self.permute = permute
        self.command = [f"gatk-{self.version}/gatk"]
        """        
    latest_version = "4.0.3"


    def fetch_version(self, version: str, target_filename: Path) -> None:
        """
        Takes care of the tool download.

        Overrides the ExternalAlgorithm methood. Downloads the external method 
        to the prebuild location specified by the corresponding
        ExternalAlgorithmStore and packs it into a tar.gz file.

        Parameters
        ----------
        version : str
            The tool version to be used.
        target_filename : Path
            The path to the local tar.gz file.
        """
        major = version[:version.rfind(".")]
        url = f"https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/gsea/software/desktop/{major}/GSEA_Linux_{version}.zip"
        print(url)
        download_file(url, open("test.zip", "wb"))
        download_zip_and_turn_into_tar_gzip(
            url, target_filename #  , chmod_x_files=[f"gatk-{version}/gatk"]
        )

    @property
    def name(self) -> str:
        """
        Returns the name of the external method for version handling.

        Overrides the ExternalAlgorithm method.

        Returns
        -------
        str
            Name of the external method.
        """
        return "GSEA"

    def build_cmd(
        self, output_directory: Optional[Path], ncores: int, arguments: List[str]
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
        return self.command + arguments

    def get_latest_version(self):
        """Returns the latest_version attribute."""
        return self.latest_version


'''

    def __repr__(self) -> str:
        return f"GSEA({self.tool}, {self.options})"

    def __str__(self) -> str:
        return f"GATK(tool = {self.tool}, options = {self.options})"

    def run(self) -> Callable:
        """
        Returns a callable that will run the pathway analysis.

        Returns
        -------
        Callable
            A callable that runs the actual analysis.
        """

        def call(
            input_files: List[List[AlignedSample]], output_file: Path, *args, **kwargs
        ):
            cmd_call = [self.caller_type]
            cmd_call.extend([str(posixpath) for posixpath in input_files])
            cmd_call.append(str(output_file))
            if self.caller_type == "somatic" or self.caller_type == "copynumber":
                if len(input_files) != 2:
                    raise ValueError(
                        f"Selected analysis {self.caller_type} needs a list of matched pileup files, got {input_files, type(input_files)}"
                    )
                cmd_call.extend(
                    [
                        "--output-snp",
                        str(output_file.parent / (output_file.stem + ".snp")),
                    ]
                )
                cmd_call.extend(
                    [
                        "--output-indel",
                        str(output_file.parent / (output_file.stem + ".indel")),
                    ]
                )
            cmd_call.extend(self.optionhandler.options_as_list_str(self.options))
            cmd_call = self.build_cmd(
                output_directory=output_file.parent,
                ncores=self.get_cores_needed(),
                arguments=cmd_call,
            )
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with Path(str(output_file) + ".varscan.log").open(
                "w"
            ) as stderr, output_file.open("w") as outp:
                stderr.write(" ".join(cmd_call) + "\n")
                try:
                    p2 = subprocess.Popen(cmd_call, stdout=outp, stderr=stderr)
                    p2.communicate()
                except subprocess.CalledProcessError:
                    print(f"Call failed : {' '.join(cmd_call)}")
                    raise
            if self.caller_type in ["somatic"]:
                with Path(output_file).open("w") as outp:
                    outp.write(
                        f'Results in \n{str(output_file.parent / (output_file.stem + ".snp"))}\n{str(output_file.parent / (output_file.stem + ".indel"))}'
                    )

        if self.preprocessor is not None:
            modifier = self.preprocessor.run_modifier()
            call = modifier(call)
        return call



        def call_gsea(self):
        """execute actual GSEA java program"""
        gsea_cmd = str(gsea_path.resolve())
        sentinel_path = self.output_path / "sentinel.txt"
        #os.chdir(str(self.output_path))
        def call():
            gsea_call = [
                "java",
                # '-Ddebug=true',
                "-cp",
                gsea_cmd,
                "-Xmx%im" % self.memory_in_mb,
                "xtools.gsea.Gsea",
                "-res",
                str(self.gct_file.name),
                "-cls",
                str(self.cls_file.name),
                "-chip",
                str(self.chp_file.name),
                # '-gmx', 'gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.all.v2.5.symbols.gmt'
                "-collapse",
                "true" if self.do_collapse else "false",
                "-mode",
                "Max_probe",
                "-norm",
                "meandiv",
                "-nperm",
                "%i" % self.n_permutation,
                "-permute",
                self.permute,
                "-rnd_type",
                "no_balance",
                "-scoring_scheme",
                "weighted",
                "-rpt_label",
                "my_analysis",
                "-metric",
                self.metric,
                "-sort",
                "real",
                "-order",
                "descending",
                "-include_only_symbols",
                "true",
                "-make_sets",
                "true",
                "-gmx",
                str(self.gmt_file.name),
                "-median",
                "false",
                "-num",
                "100",
                "-plot_top_x",
                "200",
                "-rnd_seed",
                "timestamp",
                "-save_rnd_lists",
                "false",
                "-set_max",
                "500",
                "-set_min",
                "5",
                "-zip_report",
                "false",
                "-out",
                ".",  # cwd is set to current dir...
                "-gui",
                "false",
            ]
            handle = open(sentinel_path, "wb")
            handle.write(" ".join(gsea_call).encode())
            print(" ".join(gsea_call))
            proc = subprocess.Popen(
                gsea_call,
                stdout=handle,
                stderr=subprocess.PIPE,
                shell=False,
                cwd=str(self.output_path),
            )
            stdout, stderr = proc.communicate()
            if proc.returncode != 0:
                print(stdout)
                print("stderr")
                print(stderr)

                raise ValueError(
                    "GSEA CALCULATION returncode != 0: %i\n%s\n%s"
                    % (proc.returncode, stdout, stderr)
                )

            handle.write(b"DONE")
            handle.close()

        job = (
            ppg.FileGeneratingJob(sentinel_path, call)
            .depends_on(self.dump_files())
            .depends_on(
                ppg.ParameterInvariant(
                    sentinel_path,
                    (
                        self.get_version(),
                        self.metric,
                        self.n_permutation,
                        self.do_collapse,
                        self.permute,
                    ),
                )
            )
        )
        job.cores_needed = 8
        return job



    def print_help(self) -> None:
        """
        Prints a help string for the external method.

        Calls the OptionHandler function to print a help string of the
        corresponding OPtionHandler object.
        """
        self.optionhandler.print_help()

    def print_tools(self) -> None:
        """
        Prints a list of accepted GATK tools.

        Calls the GATK tool to get a list of accepted tool commands and prints 
        them.
        
        [extended_summary]
        """
        self.run_command(["--list"])

    @property
    def multi_core(self) -> bool:
        """
        Returns True if the GATK call can use multiple cores.

        Overrides the ExternalAlgorithm method.
        """
        return "-nct" in self.options or "-nt" in self.options

    def build_cmd(
        self, output_directory: Optional[Path], ncores: int, arguments: List[str]
    ):
        """
        Returns a command as a list of strings to be passed to subprocess.

        Constructs an executable command to be passed to subprocess from the
        output directory, the number of cores to be used and a list of command
        options. Overrides the ExternalAlgorithm method and ignores the output
        directory and ncores parameters. Both are handled by the option list.

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
        if self.tool is None:
            return [str(self.path / self.command[0])] + arguments
        return [str(self.path / self.command[0]), self.tool] + arguments

    def get_cores_needed(self) -> int:
        """
        Returns the number of cores the GATK call will need.

        Overrides the ExternalAlgorithm abstract method.

        Returns
        -------
        int
            [description]
        """
        if "-nt" in self.options:
            return self.options["-nct"]
        elif "-nct" in self.options:
            return self.options["-nct"]
        return 1
'''