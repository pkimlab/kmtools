# -*- coding: utf-8 -*-
import sys

import numpy as np
import pandas as pd
import igraph


def ll2df(lofl, index, name, index_to_keep=None):
    """Convert a list of lists to a DataFrame.
    """
    df = pd.DataFrame(lofl, index=index.copy(), columns=index.copy())
    if index_to_keep is not None:
        index_to_drop = {c for c in index if c not in index_to_keep}
        df = df.drop(pd.Index(index_to_drop))
        df = df.drop(pd.Index(index_to_drop), axis=1)
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

    if ids_to_keep is not None:
        ids_to_keep = set(ids_to_keep)
        index_to_keep = set(vertex_df[vertex_df.index.isin(ids_to_keep)]['graph_idx'])
    else:
        index_to_keep = None

    edge_df['graph_idx_1'] = edge_df['id_1'].map(vertex_df['graph_idx'])
    edge_df['graph_idx_2'] = edge_df['id_2'].map(vertex_df['graph_idx'])

    g = igraph.Graph(len(vertex_df))
    g.add_edges(edge_df[['graph_idx_1', 'graph_idx_2']].apply(tuple, axis=1).tolist())

    if 'weight' in edge_df:
        print('Using weights...')
        sys.stdout.flush()
        g.es['weight'] = edge_df['weight'].values
        assert g.is_weighted()

    # === Vertex properties ===
    # degree
    print('degree')
    sys.stdout.flush()
    vertex_df['degree'] = g.degree()
    # closeness
    print('closeness')
    sys.stdout.flush()
    vertex_df['closeness'] = g.closeness()
    # betweenness
    print('betweenness')
    sys.stdout.flush()
    vertex_df['betweenness'] = g.betweenness()
    # clustering_coef
    print('clustering_coef')
    sys.stdout.flush()
    vertex_df['clustering_coef'] = g.transitivity_local_undirected()

    if cell_ppi:
        vertex_df.loc[vertex_df['closeness'] < 0.0001, 'closeness'] = np.nan
        vertex_df.loc[vertex_df['betweenness'] == 0, 'betweenness'] = np.nan
        vertex_df['betweenness'] = np.log10(vertex_df['betweenness'].values)
    return vertex_df

    # === Edge properties ===
    # edge_betweenness
    print('edge_betweenness')
    sys.stdout.flush()
    edge_df['edge_betweenness'] = g.edge_betweenness(weights='weight')
    # shortest_paths
    print('distance_min')
    sys.stdout.flush()
    lofl = g.shortest_paths(weights='weight')
    all_edge_df = ll2df(lofl, vertex_df.index, 'distance_min', index_to_keep)
    # edge_df = edge_df.merge(all_edge_df, on=['id_1', 'id_2'], how='outer')
    # similarity_jaccard
    print('similarity_jaccard')
    sys.stdout.flush()
    lofl = g.similarity_jaccard()
    df = ll2df(lofl, vertex_df.index, 'similarity_jaccard', index_to_keep)
    all_edge_df = all_edge_df.merge(df, on=['id_1', 'id_2'], how='outer')
    # similarity_inverse_log_weighted
    print('similarity_inverse_log_weighted')
    sys.stdout.flush()
    lofl = g.similarity_inverse_log_weighted()
    df = ll2df(lofl, vertex_df.index, 'similarity_inverse_log_weighted', index_to_keep)
    all_edge_df = all_edge_df.merge(df, on=['id_1', 'id_2'], how='outer')

    if cell_ppi:
        edge_df.loc[edge_df['edge_betweenness'] == 0, 'edge_betweenness'] = np.nan
        edge_df['edge_betweenness'] = np.log10(edge_df['edge_betweenness'])
    return vertex_df, edge_df, all_edge_df


if __name__ == '__main__':
    main()
