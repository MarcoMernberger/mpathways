#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""databases.py: Contains gene set databases for pathway analysis."""

from pathlib import Path
from typing import Optional, Callable, List, Dict, Tuple, Any
import pandas as pd
import pypipegraph as ppg

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def pathway_dataset(name):
    # select the correct dataset by name
    pass

# ExternalDataBase deals with versioning and downloads


class GMTDataset(ExternalDataBase):

    def __init__(self, name, version):
        pass

    def dump_gmt(self):
        # generate a GSEA usable input.gmt in a central location

    def get_gmt(self, outputfile: Path) -> ppg.FileGeneratingJob:
        pass        

