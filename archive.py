#  Copyright (c) 2020 , AnalyseSTRING
#  Author: Rohit Suratekar
#  Organisation: IIMCB, Warsaw
#  Contact: rsuratekar [at] iimcb [dot] gov [dot] pl
#
# All archived functions which are used in any publications or blog

from SecretColors import Palette

from functions import generate_network
from models import NetworkGenerator
from visualize import plot_network
from wrappers import ConfigWrapper, convert_data


def find_gene_networks():
    p = Palette()
    cw = ConfigWrapper("zebrafish")
    links2 = convert_data(cw.links)
    info2 = convert_data(cw.info)

    # remove prefix for gene IDs which don't have name assigned

    renames = {"ENSDARG00000": ""}

    # Color Mapping
    mapping = {
        "chl1a": p.blue(),
        "ptprga": p.red(shade=40),
        "ca16b": p.red(shade=40)
    }

    gene = "cdc25"  # Change here other genes of interest
    depth = 1  # Change here when looking for multi-level connections
    g = generate_network(gene,
                         info=info2,
                         links=links2,
                         threshold=0.4,
                         depth=depth)
    plot_network(g,
                 replace_names=renames,
                 color_mapping=mapping)


def other_networks():
    p = Palette()
    cw = ConfigWrapper("mouse")  # Change human or mouse
    links2 = convert_data(cw.links)
    info2 = convert_data(cw.info)

    # remove prefix for gene IDs which don't have name assigned

    renames = {"ENSDARG00000": ""}

    # Color Mapping
    mapping = {
        "Chl1": p.blue(),
        "Ptprg": p.red(shade=40),
        "CA16": p.red(shade=40)
    }

    gene = "Chl1"  # Change here other genes of interest
    depth = 1  # Change here when looking for multi-level connections
    g = generate_network(gene,
                         info=info2,
                         links=links2,
                         threshold=0.4,
                         depth=depth)
    plot_network(g,
                 replace_names=renames,
                 color_mapping=mapping)


def calculate_score_statistics():
    cw = ConfigWrapper("mouse")  # Change human or mouse
    links2 = convert_data(cw.links)
    info2 = convert_data(cw.info)

    gene1 = "Chl1"  # Change here other genes of interest
    gene2 = "Ptprg"  # Change here other genes of interest

    n = NetworkGenerator(info=info2, links=links2)
    n.show_entry(gene1, gene2)


def run():
    calculate_score_statistics()
