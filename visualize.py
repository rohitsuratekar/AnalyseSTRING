#  Copyright (c) 2020 , AnalyseSTRING
#  Author: Rohit Suratekar
#  Organisation: IIMCB, Warsaw
#  Contact: rsuratekar [at] iimcb [dot] gov [dot] pl
#
# All visualization functions

import numpy as np
from SecretColors import Palette
from SecretPlots import NetworkPlot


def paths_to_connections(data, paths):
    connections = []
    for path in paths:
        for i in range(1, len(path)):
            start = path[i - 1]
            end = path[i]
            row = np.array([start, end])
            ind = np.where(np.all(data == row, axis=1))
            if len(data[ind]) == 0:
                connections.append([end, start, 1])
            else:
                connections.append([start, end, 1])

    return connections


def plot_network(data, max_cols=None, show=True, color_mapping: dict = None):
    p = Palette()
    network_data = []
    if len(data) > len(p.get_color_list):
        colors = p.random(no_of_colors=len(data))
    else:
        colors = p.get_color_list[len(data)]
        if len(data) == 1:
            colors = [colors]

    mapping_colors = {}

    for pair in data:
        mapping_colors[pair[0]] = colors[len(mapping_colors)]
        for gene in pair[1]:
            network_data.append([pair[0], gene, 1])

    n = NetworkPlot(network_data)
    if max_cols is not None:
        n.max_columns = max_cols
    n.colors = p.gray(shade=30)
    if color_mapping is None:
        n.colors_mapping = mapping_colors
    else:
        n.colors_mapping = color_mapping
    n.line_decoration = False
    if show:
        n.show()
    return n


def plot_connections(data, max_cols=None, show=True,
                     color_mapping: dict = None):
    p = Palette()
    all_nodes = []
    for d in data:
        for i, node in enumerate(d):
            if i == 0:
                all_nodes.append([node, d[i + 1], 1])
            elif i == len(d) - 1 and len(d) == 2:
                continue
            else:
                all_nodes.append([d[i - 1], node, 1])

    mapping_colors = {
        data[0][0]: p.blue(),
        data[0][-1]: p.red()
    }

    n = NetworkPlot(all_nodes)
    if max_cols is not None:
        n.max_columns = max_cols
    n.colors = p.gray(shade=30)
    if color_mapping is not None:
        n.colors_mapping = color_mapping
    else:
        n.colors_mapping = mapping_colors
    if show:
        n.show()
    return n
