#  Copyright (c) 2020 , AnalyseSTRING
#  Author: Rohit Suratekar
#  Organisation: IIMCB, Warsaw
#  Contact: rsuratekar [at] iimcb [dot] gov [dot] pl
#
#  All helper wrappers will go in this file

import json

import pandas as pd
from pandas.errors import ParserError


class ConfigWrapper:
    """
    Simple config file wrapper which can be used to automatically extract
    various filenames from the config file.
    """
    ORGANISMS = ["zebrafish", "human", "rat", "mouse"]

    def __init__(self, organism: str, filename: str = None):
        if filename is None:
            filename = "config.json"

        with open(filename) as f:
            self.data = json.load(f)

        self._organism = organism

    @property
    def organism(self) -> str:
        return self._organism.strip().lower()

    def _validate_organism(self):
        if self.organism not in ConfigWrapper.ORGANISMS:
            raise ValueError(f"Invalid organism. Currently supported "
                             f"organisms: {ConfigWrapper.ORGANISMS}")

    @property
    def actions(self) -> str:
        self._validate_organism()
        return self.data[self.organism]["actions"]

    @property
    def links(self) -> str:
        self._validate_organism()
        return self.data[self.organism]["full_links"]

    @property
    def info(self) -> str:
        self._validate_organism()
        return self.data[self.organism]["info"]


def convert_data(filename: str) -> pd.DataFrame:
    """
    Converts STRING database files into pandas DataFrame.
    This assumes, STRING database files are tab separated.
    :param filename: Full path of the file
    :return: Pandas DataFrame
    """
    try:
        return pd.read_table(filename, delim_whitespace=True)
    except ParserError:
        # Info files are tab delimited
        return pd.read_csv(filename, delimiter="\t")


def run():
    cw = ConfigWrapper("zebrafish")
    print(cw.info)
