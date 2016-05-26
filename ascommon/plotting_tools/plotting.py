import scipy as sp
import seaborn as sns


def make_plot_with_corr(x, y, data, ax, corr_type='spearman'):
    sns.regplot(
        x=x, y=y, data=data,
        scatter_kws={
            'alpha': 0.8,
            'linewidth': 1,
            'edgecolor': 'w',
            's': 50},
        line_kws={'color': 'k', 'linewidth': 2},
        ax=ax
    )
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xpos = xlim[0] + 0.02 * (xlim[1] - xlim[0])
    ypos = ylim[0] + 0.98 * (ylim[1] - ylim[0])
    # Show correlation
    if corr_type == 'spearman':
        corr, pvalue = sp.stats.spearmanr(data[[x, y]].dropna()[x], data[[x, y]].dropna()[y])
    elif corr_type == 'pearson':
        corr, pvalue = sp.stats.pearsonr(data[[x, y]].dropna()[x], data[[x, y]].dropna()[y])
    ax.text(
        xpos, ypos, 'R: {:.3f}\np-value: {:.1e}'.format(corr, pvalue), fontsize=18,
        horizontalalignment='left',
        verticalalignment='top',
    )
