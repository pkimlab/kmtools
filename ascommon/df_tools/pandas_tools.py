import numpy as np


def split_df(df, n_chunks):
    """Split DataFrame `df` into `nchunks` chunks."""
    chunk_size = int(np.ceil(df.shape[0] / n_chunks))
    assert n_chunks * chunk_size >= df.shape[0]
    chunks = []
    for i in range(0, df.shape[0], chunk_size):
        chunks.append(df[i:i + chunk_size])
    assert len(chunks) == n_chunks
    return chunks
