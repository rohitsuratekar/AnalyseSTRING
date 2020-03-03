#  Copyright (c) 2020 , AnalyseSTRING
#  Author: Rohit Suratekar
#  Organisation: IIMCB, Warsaw
#  Contact: rsuratekar [at] iimcb [dot] gov [dot] pl
#
#  All object models and classes

import pandas as pd

from constants import *


class Score:
    """
    Simple class to hold the column information
    """

    def __init__(self, name):
        self.name = name
        self.use = False


class NetworkFinder:
    """
    Class used to calculate the network from STRING database
    """

    def __init__(self, *, info, links):
        self._info = None
        self._links = None
        self._raw_info = info
        self._raw_links = links
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

        self.experiments.use = True
        self.experiments_trans.use = True

        self.threshold = 0.4

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

    @property
    def ids(self):
        """
        :return: Dictionary or Mapping of Protein IDs (as index) and
        preferred name (as value)
        """
        if self._info is None:
            k = self._raw_info
            k = pd.Series(k[PREFERRED_NAME].values, index=k[
                PROTEIN_EXTERNAL_ID].values)
            self._info = k.to_dict()
        return self._info

    @property
    def interactions(self) -> pd.DataFrame:
        """
        Protein-Protein interactions and their individual scores.
        """
        if self._links is None:
            # Get data
            self._links = self._raw_links
            # Convert IDs to regular names
            self._links[PROTEIN_1] = self._links[PROTEIN_1].map(
                lambda x: self.ids[x])
            self._links[PROTEIN_2] = self._links[PROTEIN_2].map(
                lambda x: self.ids[x])

            # Calculate temporary score for given traits
            temp = "temp"
            self._links[temp] = 0
            for c in self.columns:
                if c.use:
                    self._links[temp] = self._links[temp] + self._links[c.name]

            # Take only entries which have more score than the user set
            # threshold (database values are multiplied by 1000)
            self._links = (self._links[
                self._links[temp] >= self.threshold * 1000].reset_index(
                drop=True))

            # Remove the temporary column
            del self._links[temp]

        return self._links

    def find_interactors(self, gene: str):
        """
        Finds interacting proteins with given settings
        :param gene: Gene name
        :return: List of interacting partners
        """
        data = self.interactions[
            self.interactions[PROTEIN_1].str.lower() == gene.strip(
            ).lower()]
        return list(data[PROTEIN_2].values)

    def generate_network(self, genes, *, level: int = 1):
        # Sanity check
        if level < 1 or type(level) != int:
            raise Exception(
                "Invalid Level. Please use integers greater than 0")
        # if only single gene is given, convert it into a list
        if type(genes) == str:
            genes = [genes]
        net_stat = []
        all_genes = []
        secondary_genes = []
        # Start network expansion
        for i in range(level):
            temp = []
            if len(net_stat) == 0:
                secondary_genes = genes
            for g in secondary_genes:
                if g not in all_genes:
                    links = self.find_interactors(g)
                    net_stat.append((g, links))
                    all_genes.append(g)
                    temp.extend(links)
            secondary_genes = temp
            all_genes = list(set(all_genes))

        return net_stat
