#  Copyright (c) 2020 , AnalyseSTRING
#  Author: Rohit Suratekar
#  Organisation: IIMCB, Warsaw
#  Contact: rsuratekar [at] iimcb [dot] gov [dot] pl
#
# Simple scripts to analyse the STRING database (https://string-db.org/)


from collections import defaultdict

import pandas as pd
from SecretColors import Palette
from SecretPlots import NetworkPlot
from pandas.errors import ParserError

from models import NetworkFinder
from wrappers import ConfigWrapper


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


def convert_to_digital(network):
    digital = []
    for gene_tuple in network:
        for g in gene_tuple[1]:
            g2 = gene_tuple[0]
            # Remove the gene ID
            if "ENSDARG" in str(g):
                g = g.replace("ENSDARG", "").replace("0", "")
            if "ENSDARG" in str(g2):
                g2 = g2.replace("ENSDARG", "").replace("0", "")
            digital.append([g2, g, 1])
    return digital


def find_connection(gene1, gene2, info, links):
    n = NetworkFinder(info=info, links=links)
    if len(n.find_interactors(gene1)) == 0:
        print(f"{gene1} has no interactors in current database")
        return
    if len(n.find_interactors(gene2)) == 0:
        print(f"{gene2} has no interactors in current database")
        return
    all_genes = list(n.ids.values())
    back_track = defaultdict(list)
    done_genes = []
    next_genes = []
    while len(all_genes) != 0:
        if len(next_genes) == 0:
            next_genes = [gene1]
        else:
            next_genes = []

        tmp_list = []
        found_gene = False
        for gene in next_genes:
            if gene not in done_genes:
                k = n.find_interactors(gene)
                if gene2 in k:
                    found_gene = True
                for t in k:
                    back_track[t].append(gene)

                tmp_list.extend(k)
                done_genes.append(gene)
                all_genes.remove(gene)

        next_genes = list(set(tmp_list))
        if found_gene:
            break
        print(next_genes)

    print(back_track)


def visualize_network(gene, *, info, links, level=1):
    p = Palette()
    n = NetworkFinder(info=info, links=links)
    k = n.generate_network(gene, level=level)
    data = convert_to_digital(k)
    np = NetworkPlot(data)
    np.node_gap = 2
    np.node_width = 3
    # np.add_text_options(fontsize=9)
    np.colors = p.gray(shade=30)
    np.colors_mapping = {"ptprz1a": p.blue(), "chl1a": p.green(),
                         "ptprga": p.blue(), "ca16b": p.blue(),
                         "cntn1a": p.blue(), gene: p.red()}
    # np.save("network.png", dpi=300, format="png")
    np.show()


def run():
    cw = ConfigWrapper("zebrafish")
    links = convert_data(cw.links)
    info = convert_data(cw.info)

    # Find connections between two genes
    # This will take long time if connections are distinctly connected and
    # number of gene available
    # BE SURE that gene provided are present in your database else Divisible
    # by Zero error will be raised by NetWork visualizer from SecretPlot
    # library
    find_connection("nkx2.5", "gata4", info=info, links=links)


if __name__ == "__main__":
    run()
