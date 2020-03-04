#  Copyright (c) 2020 , AnalyseSTRING
#  Author: Rohit Suratekar
#  Organisation: IIMCB, Warsaw
#  Contact: rsuratekar [at] iimcb [dot] gov [dot] pl
#
# Simple scripts to analyse the STRING database (https://string-db.org/)
#
# Important Note about the dataset files: Use dataset files from the same
# release i.e. download ALL the files from same date.

import pprint
from typing import Union

import pandas as pd

from constants import PROTEIN_1, PROTEIN_2
from models import NetworkGenerator, ClusterNode
from visualize import plot_connections
from wrappers import ConfigWrapper, convert_data


def _sanity_check(threshold):
    if threshold > 1 or threshold < 0:
        raise ValueError("Please use threshold between 0 to 1")


def generate_network(genes: Union[str, list], *,
                     depth: int = 1,
                     threshold: float = 0.4,
                     info: pd.DataFrame,
                     links: pd.DataFrame,
                     network: NetworkGenerator = None,
                     print_network: bool = True) -> list:
    # Sanity check
    _sanity_check(threshold)
    if network is None:
        network = NetworkGenerator(info=info, links=links)
        network.experiments.use = True
        network.experiments_trans.use = True

    network.threshold = threshold

    nr = network.generate_network(genes, depth=depth)

    if print_network:
        pprint.pprint(nr)

    return nr


def generate_cluster(data):
    clusters = []
    for d in data:
        values = [len(set(x).intersection(d)) for x in clusters]
        if sum(values) == 0:
            clusters.append(list(d))
        else:
            tmp = list(d)
            tmp2 = []
            for i in range(len(values)):
                if values[i] != 0:
                    tmp.extend(clusters[i])
                else:
                    tmp2.append(clusters[i])

            tmp2.append(list(set(tmp)))
            clusters = tmp2

    return clusters


def _calculate_distance(data, start_node):
    all_nodes = list(set(data.flatten()))
    all_nodes = {x: ClusterNode(x) for x in all_nodes}
    all_nodes[start_node].value = 0
    counter = 1
    while counter != len(all_nodes):
        temp_count = 0
        for d in data:
            d1 = all_nodes[d[0]]  # type:ClusterNode
            d2 = all_nodes[d[1]]  # type:ClusterNode
            if d2.add_distance(d1.value + 1, d1.name):
                temp_count += 1
            if d1.add_distance(d2.value + 1, d2.name):
                temp_count += 1

        if temp_count == 0:
            break
        else:
            counter += temp_count

    return all_nodes


def find_shortest_connection(gene1, gene2, *,
                             links: pd.DataFrame,
                             info: pd.DataFrame,
                             network: NetworkGenerator = None,
                             threshold: float = 0.4,
                             print_paths: bool = True):
    if network is None:
        data = NetworkGenerator(info=info, links=links)
        data.experiments.use = True
        data.experiments_trans.use = True
    else:
        data = network

    _sanity_check(threshold)
    data.threshold = threshold
    names = {data.names[x]: x for x in data.names}
    data = data.get_filtered()
    data = data[[PROTEIN_1, PROTEIN_2]].applymap(lambda x: names[x]).values

    if gene2 not in data.flatten():
        raise Exception(f"{gene2} is not present in the database. Try "
                        f"reducing threshold value?")
    d = _calculate_distance(data, gene1)
    all_paths = []

    def _find_path(data_list, temp_ng):
        temp_path = data_list
        while True:
            extra_paths = []
            if temp_ng is None:
                break
            ng = d[temp_ng]  # type: ClusterNode
            if ng is not None:
                if len(ng.alternatives) > 0:
                    alt = [x[1] for x in ng.alternatives]
                    extra_paths.extend(alt)
                    temp_path.append(ng.name)
                else:
                    temp_path.append(ng.name)
                temp_ng = ng.closest
            else:
                temp_path.append(ng.name)

            for e in extra_paths:
                ex_path = [x for x in temp_path]
                _find_path(ex_path, e)

        temp_path = list(reversed(temp_path))
        all_paths.append(temp_path)

    _find_path([], gene2)
    if print_paths:
        for p in all_paths:
            current_str = ""
            for k in p:
                current_str += f"{k} -> "
            pprint.pprint(current_str[:-4])

    if len(all_paths) == 1:
        if len(all_paths[0]) == 1:
            raise Exception(
                f"Unfortunately, there is no possible connection "
                f"from '{gene1}' to '{gene2}' with given threshold "
                f"'{threshold}'. Try reducing threshold.")

    return all_paths


if __name__ == "__main__":
    cw = ConfigWrapper("zebrafish")
    links2 = convert_data(cw.links)
    info2 = convert_data(cw.info)

    g = find_shortest_connection("nkx2.5", "mef2ca", info=info2, links=links2)
    # g = generate_network("slc35a5", links=links2, info=info2, threshold=0.2,
    #                      print_network=False)
    # plot_network(g)
    plot_connections(g)
