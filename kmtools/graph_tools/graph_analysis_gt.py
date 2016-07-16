# -*- coding: utf-8 -*-
import sys

import numpy as np
import pandas as pd
# import igraph
import graph_tool.all as gt


def ll2df(lofl, index, name):
    """Convert a list of lists to a DataFrame.
    """
    df = pd.DataFrame(lofl, index=index.copy(), columns=index.copy())
    df.index.names = ['id_2']
    df.columns.names = ['id_1']
    df = df.unstack().reset_index().rename(columns={0: name})
    df.loc[np.isinf(df[name]), name] = np.nan
    return df


def main(edge_df, cell_ppi=False, ids_to_keep=None):
    """Compute vertex and edge properties for a graph defined by `edge_df`.
    """
    vertex_df = (
        pd.concat([edge_df['id_1'], edge_df['id_2']], ignore_index=True)
        .drop_duplicates()
        .reset_index()
        .rename(columns={'index': 'graph_idx', 0: 'id'})
        .set_index('id')
    )
    vertex_df['graph_idx'] = range(len(vertex_df))

    edge_df['graph_idx_1'] = edge_df['id_1'].map(vertex_df['graph_idx'])
    edge_df['graph_idx_2'] = edge_df['id_2'].map(vertex_df['graph_idx'])

    # g = igraph.Graph(len(vertex_df))
    # g.add_edges(edge_df[['graph_idx_1', 'graph_idx_2']].apply(tuple, axis=1).tolist())

    g = gt.Graph(directed=False)
    g.add_vertex(len(vertex_df))
    g.add_edge_list(edge_df[['graph_idx_1', 'graph_idx_2']].values)
    # for a, b in edge_df[['graph_idx_1', 'graph_idx_2']].values:
    #     g.add_edge(g.vertex(a), g.vertex(b))

    if 'weight' in edge_df:
        print('Using weights...')
        sys.stdout.flush()
        g_weight = g.new_edge_property('double', vals=edge_df['weight'].values)
    else:
        g_weight = None

#    # === Vertex properties ===
#    # degree
#    print('degree')
#    sys.stdout.flush()
#    vertex_df['degree'] = g.degree_property_map(deg='total', weight=g_weight).get_array()
#    # closeness
#    print('closeness')
#    sys.stdout.flush()
#    vertex_df['closeness'] = gt.closeness(g, weight=g_weight).get_array()
#    # Pagerank
#    print('pagerank')
#    sys.stdout.flush()
#    vertex_df['pagerank'] = gt.pagerank(g, weight=g_weight).get_array()
#    # clustering_coef
#    print('clustering_coef')
#    sys.stdout.flush()
#    vertex_df['clustering_coef'] = gt.local_clustering(g).get_array()
#
#    if cell_ppi:
#        vertex_df.loc[vertex_df['closeness'] < 0.0001, 'closeness'] = np.nan
#        vertex_df.loc[vertex_df['betweenness'] == 0, 'betweenness'] = np.nan
#        vertex_df['betweenness'] = np.log10(vertex_df['betweenness'].values)
#
#    # === Edge properties ===
#    # betweenness
#    print('betweenness')
#    sys.stdout.flush()
#    betweenness = gt.betweenness(g, weight=g_weight)
#    vertex_df['betweenness'] = betweenness[0].get_array()
#    edge_df['edge_betweenness'] = betweenness[1].get_array()

    # === All by all ===
    if ids_to_keep is not None:
        ids_to_keep = set(ids_to_keep)
        index_to_keep = list(set(vertex_df[vertex_df.index.isin(ids_to_keep)]['graph_idx']))
    else:
        index_to_keep = list(vertex_df.index)

    # shortest_paths
    print('distance_min')
    sys.stdout.flush()
    lofl = (
        gt
        .shortest_distance(g, weights=g_weight)
        .get_2d_array(index_to_keep)
        [:, index_to_keep]
    )
    df = ll2df(lofl, pd.Index(index_to_keep), 'distance_min')
    all_edge_df = df

    # trust_transitivity
    print('trust_transitivity')
    sys.stdout.flush()
    lofl = (
        gt
        .trust_transitivity(g, trust_map=g_weight)
        .get_2d_array(index_to_keep)
        [:, index_to_keep]
    )
    df = ll2df(lofl, pd.Index(index_to_keep), 'trust_transitivity')
    all_edge_df = all_edge_df.merge(df, on=['id_1', 'id_2'], how='outer')

    return vertex_df, edge_df, all_edge_df

#    # edge_df = edge_df.merge(all_edge_df, on=['id_1', 'id_2'], how='outer')
#    # similarity_jaccard
#    print('similarity_jaccard')
#    sys.stdout.flush()
#    lofl = g.similarity_jaccard()
#    df = ll2df(lofl, vertex_df.index, 'similarity_jaccard', index_to_keep)
#    all_edge_df = all_edge_df.merge(df, on=['id_1', 'id_2'], how='outer')
#
#    # similarity_inverse_log_weighted
#    print('similarity_inverse_log_weighted')
#    sys.stdout.flush()
#    lofl = g.similarity_inverse_log_weighted()
#    df = ll2df(lofl, vertex_df.index, 'similarity_inverse_log_weighted', index_to_keep)
#    all_edge_df = all_edge_df.merge(df, on=['id_1', 'id_2'], how='outer')
#
#    if cell_ppi:
#        edge_df.loc[edge_df['edge_betweenness'] == 0, 'edge_betweenness'] = np.nan
#        edge_df['edge_betweenness'] = np.log10(edge_df['edge_betweenness'])
#
#    return vertex_df, edge_df, all_edge_df


if __name__ == '__main__':
    main()
