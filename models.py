#  Copyright (c) 2020 , AnalyseSTRING
#  Author: Rohit Suratekar
#  Organisation: IIMCB, Warsaw
#  Contact: rsuratekar [at] iimcb [dot] gov [dot] pl
#
#  All object models and classes

from typing import Union

import numpy as np

from constants import *
from wrappers import *


class Score:
    """
    Simple class to hold the column information
    """

    def __init__(self, name):
        self.name = name
        self.use = False


class NetworkGenerator:
    def __init__(self, *, info: pd.DataFrame, links: pd.DataFrame):
        self._info = info
        self._links = links
        self._names = None

        self.neighborhood = Score(NEIGHBORHOOD)
        self.neighborhood_trans = Score(NEIGHBORHOOD_TRANSFERRED)
        self.fusion = Score(FUSION)
        self.cooccurence = Score(COOCCURENCE)
        self.homology = Score(HOMOLOGY)
        self.coexpression = Score(COEXPRESSION)
        self.coexpression_trans = Score(COEXPRESSION_TRANSFERRED)
        self.experiments = Score(EXPERIMENTS)
        self.experiments_trans = Score(EXPERIMENTS_TRANSFERRED)
        self.text_mining = Score(TEXT_MINING)
        self.text_mining_trans = Score(TEXT_MINING_TRANSFERRED)
        self.threshold = 0.4

    @property
    def info(self) -> pd.DataFrame:
        return self._info

    @property
    def links(self) -> pd.DataFrame:
        return self._links

    @property
    def names(self) -> dict:
        if self._names is None:
            self._names = self._attach_names()
        return self._names

    @property
    def columns(self):
        """
        :return: List of all columns
        """
        return [self.neighborhood, self.neighborhood_trans, self.fusion,
                self.cooccurence, self.homology, self.coexpression,
                self.coexpression_trans, self.experiments,
                self.experiments_trans, self.text_mining,
                self.text_mining_trans]

    def _attach_names(self):
        data = self.info.set_index(PROTEIN_EXTERNAL_ID)
        data = data[[PREFERRED_NAME]]
        data = data.T.to_dict()
        data = {data[x][PREFERRED_NAME]: x for x in data}
        return data

    def get_filtered(self) -> pd.DataFrame:
        df = self.links

        # Calculate temporary score for given traits
        temp = "temp"
        df[temp] = 0
        for c in self.columns:
            if c.use:
                df[temp] = df[temp] + df[c.name]

        # Take only entries which have more score than the user set
        # threshold (database values are multiplied by 1000)
        df = (df[df[temp] >= self.threshold * 1000].reset_index(
            drop=True))

        del df[temp]

        return df

    @staticmethod
    def _find_connections(df, gene):
        """
          Finds interacting proteins with given settings
          :param gene: Gene name
          :return: List of interacting partners
          """
        data = df[df[PROTEIN_1] == gene]
        return list(data[PROTEIN_2].values)

    def generate_network(self, genes: Union[str, list], *, depth: int = 1):
        if type(genes) == str:
            genes = [genes]
        # Sanity checks
        for x in genes:
            if x not in self.names.keys():
                raise ValueError(f"'{genes}' is not present in the links data "
                                 f"provided. Please provide EXACT name ("
                                 f"case-sensitive) present in the STRING "
                                 f"dataset")

        if depth < 1 or type(depth) != int:
            raise Exception(
                "Invalid Level. Please use integers greater than 0")

        # Convert all genes into the their names
        genes = [self.names[x] for x in genes]
        all_links = self.get_filtered()

        net_stat = []
        all_genes = []
        secondary_genes = []
        # Start network expansion
        for i in range(depth):
            temp = []
            if len(net_stat) == 0:
                secondary_genes = genes
            for g in secondary_genes:
                if g not in all_genes:
                    links = self._find_connections(all_links, g)
                    net_stat.append((g, links))
                    all_genes.append(g)
                    temp.extend(links)
            secondary_genes = temp
            all_genes = list(set(all_genes))

        gene_mapping = {self.names[x]: x for x in self.names}
        renamed = []
        for pair in net_stat:
            temp_connections = pair[1]
            temp_connections = [gene_mapping[x] for x in temp_connections]
            renamed.append((gene_mapping[pair[0]], temp_connections))

        return renamed


class ClusterNode:
    def __init__(self, name):
        self.name = name
        self.value = np.inf
        self.closest = None
        self.alternatives = []

    def add_distance(self, value, node):
        if value < self.value:
            self.value = value
            self.closest = node
            return True
        elif value == self.value and value != np.inf:
            new_alt = []
            for a in self.alternatives:
                if a[0] < value:
                    new_alt.append(a)
            if self.closest != node:
                new_alt.append((value, node))
            self.alternatives = new_alt

        return False

    def __repr__(self):
        return f"Node( {self.value} )"


def run():
    cw = ConfigWrapper("zebrafish")
    links = convert_data(cw.links)
    info = convert_data(cw.info)

    ng = NetworkGenerator(info=info, links=links)
    ng.experiments.use = True
    ng.experiments_trans.use = True
