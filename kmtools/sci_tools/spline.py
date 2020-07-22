import numpy as np
from scipy import stats


def get_regline(x, y, bins=1_000):
    """Preprocess large volumes of data for :class:`scipy.interpolate.splrep`.

    Examples:
        >>> from scipy import interpolate

        >>> x = np.linspace(0, np.pi * 6, 10_000)
        >>> y = np.sin(x) + np.random.randn(*x.shape) * 0.1
        >>> xyw = get_regline(x, y, bins=10)

        >>> tck = interpolate.splrep(*xyw)

        >>> xnew = np.linspace(x.min(), x.max(), 100)
        >>> ynew = interpolate.splev(xnew, tck)

        >>> plt.plot(xyw[0], xyw[1], '.')
        >>> plt.plot(xnew, ynew)
    """
    mean, count = [
        stats.binned_statistic(x, y, statistic, bins=bins, range=(x.min(), x.max()))
        for statistic in ["mean", "count"]
    ]
    assert (mean.bin_edges == count.bin_edges).all()

    include_mask = ~np.isnan(mean.statistic)

    x = (mean.bin_edges[:-1][include_mask] + mean.bin_edges[1:][include_mask]) / 2
    y = mean.statistic[include_mask]
    w = count.statistic[include_mask]

    return x, y, w
