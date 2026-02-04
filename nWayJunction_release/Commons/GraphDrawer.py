#!/usr/bin/python
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import random
from pathlib import Path


def draw_graph(list_of_records):
    record = list_of_records[4]
    junction = record.list_of_junctions[0]
    list_from = []
    list_to = []
    for ranges in range(len(junction.lenghts_of_segments)):
        if ranges == 0:
            list_from.append('Stem' + str(ranges))
            list_to.append('Stem' + str(ranges + 1))
        elif ranges == len(junction.lenghts_of_segments)-1:
            list_from.append('Stem' + str(0))
            list_to.append('Stem' + str(ranges))
        else:
            list_from.append('Stem' + str(ranges))
            list_to.append('Stem' + str(ranges + 1))

    df = pd.DataFrame({ 'from':list_from, 'to':list_to}, index = [1, 2, 3, 4])

    G=nx.from_pandas_edgelist(df, 'from', 'to')

    nx.draw(G, with_labels=True, node_size=1500, node_color="skyblue", node_shape="s", alpha=0.5, linewidths=40)

    plt.show()


def draw_graph_2(list_of_records):
    record = list_of_records[4]
    junction = record.list_of_junctions[0]
    edges = []
    labels_egde = {}
    degree = u"\u00b0"
    for ranges in range(len(junction.lenghts_of_segments)):
        if ranges == 0:
            edges.append(['Stem' + str(ranges+1),'Stem' + str(ranges + 2)])
            labels_egde[('Stem' + str(ranges+1),'Stem' + str(ranges + 2))] = str(str(junction.lenghts_of_segments[ranges]) + '\n' + str(junction.list_of_segment_seq[ranges]) + '\n' + str(junction.list_of_segment_db[ranges]) + '\n' + str('('+ str(random.randrange(0,178))  +', '+ str(random.randrange(0,178)) +  ', ' + str(random.randrange(0,178)) +  ')'))
        elif ranges == len(junction.lenghts_of_segments)-1:
            edges.append(['Stem' + str(1),'Stem' + str(ranges+1)])
            labels_egde[('Stem' + str(1),'Stem' + str(ranges+1))] = str(str(junction.lenghts_of_segments[ranges]) + '\n' + str(junction.list_of_segment_seq[ranges]) + '\n' + str(junction.list_of_segment_db[ranges]) + '\n' + str('('+ str(random.randrange(0,178)) +  ', ' 
+ str(random.randrange(0,178)) +  ', ' + str(random.randrange(0,178)) +   ')'))
        else:
            edges.append(['Stem' + str(ranges+1), 'Stem' + str(ranges + 2)])
            labels_egde[('Stem' + str(ranges+1), 'Stem' + str(ranges + 2))] = str(str(junction.lenghts_of_segments[ranges]) + '\n' + str(junction.list_of_segment_seq[ranges]) + '\n' + str(junction.list_of_segment_db[ranges]) + '\n' + str('('+ str(random.randrange(0,178)) +  ', '+ str(random.randrange(0,178)) +  ', '+ str(random.randrange(0,178)) +   ')'))

    G=nx.Graph()
    G.add_edges_from(edges)
    pos = nx.spring_layout(G)
    plt.figure()

    nx.draw(G,pos,edge_color='black',width=1,linewidths=1,\
    node_size=2150,node_color='lightgrey',alpha=0.9,\
    labels={node:node for node in G.nodes()})


    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels_egde,font_color='black', font_size= 20 )
    plt.axis('off')
    plt.savefig(Path('.\output\demo.png'), transparent=True)
    plt.show()
