#! /usr/bin/env python3
# coding=utf-8

#import scanpy as sc

import collections.abc as cabc
from itertools import product
from typing import Optional, Union, Mapping  # Special
from typing import Sequence, Collection, Iterable  # ABCs
from typing import Tuple, List  # Classes

import numpy as np
import pandas as pd
from anndata import AnnData
from cycler import Cycler
from matplotlib.axes import Axes
from pandas.api.types import is_categorical_dtype
from scipy.sparse import issparse
from matplotlib import pyplot as pl
from matplotlib import rcParams
from matplotlib import gridspec
from matplotlib import patheffects
from matplotlib.colors import is_color_like, Colormap, ListedColormap

from scanpy._compat import Literal
from scanpy import get
from scanpy import logging as logg
from scanpy._settings import settings
from scanpy._utils import sanitize_anndata, _doc_params
from scanpy._compat import Literal
from scanpy.plotting import _utils
from scanpy.plotting._utils import scatter_base, scatter_group, setup_axes
from scanpy.plotting._utils import ColorLike, _FontWeight, _FontSize
from scanpy.plotting._docs import doc_scatter_basic, doc_show_save_ax, doc_common_plot_args

_VarNames = Union[str, Sequence[str]]


def check_var_names_type(var_names, var_group_labels, var_group_positions):
    """
    checks if var_names is a dict. Is this is the cases, then set the
    correct values for var_group_labels and var_group_positions

    Returns
    -------
    var_names, var_group_labels, var_group_positions

    """
    if isinstance(var_names, cabc.Mapping):
        if var_group_labels is not None or var_group_positions is not None:
            logg.warning(
                "`var_names` is a dictionary. This will reset the current "
                "value of `var_group_labels` and `var_group_positions`."
            )
        var_group_labels = []
        _var_names = []
        var_group_positions = []
        start = 0
        for label, vars_list in var_names.items():
            if isinstance(vars_list, str):
                vars_list = [vars_list]
            # use list() in case var_list is a numpy array or pandas series
            _var_names.extend(list(vars_list))
            var_group_labels.append(label)
            var_group_positions.append((start, start + len(vars_list) - 1))
            start += len(vars_list)
        var_names = _var_names

    elif isinstance(var_names, str):
        var_names = [var_names]

    return var_names, var_group_labels, var_group_positions

def prepare_dataframe(
    adata: AnnData,
    var_names: Union[_VarNames, Mapping[str, _VarNames]],
    groupby: Optional[str] = None,
    use_raw: Optional[bool] = None,
    log: bool = False,
    num_categories: int = 7,
    layer=None,
    gene_symbols: Optional[str] = None,
):
    """
    Given the anndata object, prepares a data frame in which the row index are the categories
    defined by group by and the columns correspond to var_names.

    Parameters
    ----------
    adata
        Annotated data matrix.
    var_names
        `var_names` should be a valid subset of  `adata.var_names`.
    groupby
        The key of the observation grouping to consider. It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories`.
    use_raw
        Use `raw` attribute of `adata` if present.
    log
        Use the log of the values
    num_categories
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.
    gene_symbols
        Key for field in .var that stores gene symbols.

    Returns
    -------
    Tuple of `pandas.DataFrame` and list of categories.
    """
    from scipy.sparse import issparse

    sanitize_anndata(adata)
    if use_raw is None and adata.raw is not None:
        use_raw = True
    if isinstance(var_names, str):
        var_names = [var_names]

    if groupby is not None:
        if groupby not in adata.obs_keys():
            raise ValueError(
                'groupby has to be a valid observation. '
                f'Given {groupby}, valid observations: {adata.obs_keys()}'
            )

    if gene_symbols is not None and gene_symbols in adata.var.columns:
        # translate gene_symbols to var_names
        # slow method but gives a meaningful error if no gene symbol is found:
        translated_var_names = []
        for symbol in var_names:
            if symbol not in adata.var[gene_symbols].values:
                logg.error(
                    f"Gene symbol {symbol!r} not found in given "
                    f"gene_symbols column: {gene_symbols!r}"
                )
                return
            translated_var_names.append(
                adata.var[adata.var[gene_symbols] == symbol].index[0]
            )
        symbols = var_names
        var_names = translated_var_names
    if layer is not None:
        if layer not in adata.layers.keys():
            raise KeyError(
                f'Selected layer: {layer} is not in the layers list. '
                f'The list of valid layers is: {adata.layers.keys()}'
            )
        matrix = adata[:, var_names].layers[layer]
    elif use_raw:
        matrix = adata.raw[:, var_names].X
    else:
        matrix = adata[:, var_names].X

    if issparse(matrix):
        matrix = matrix.toarray()
    if log:
        matrix = np.log1p(matrix)

    obs_tidy = pd.DataFrame(matrix, columns=var_names)
    if groupby is None:
        groupby = ''
        categorical = pd.Series(np.repeat('', len(obs_tidy))).astype('category')
    else:
        if not is_categorical_dtype(adata.obs[groupby]):
            # if the groupby column is not categorical, turn it into one
            # by subdividing into  `num_categories` categories
            categorical = pd.cut(adata.obs[groupby], num_categories)
        else:
            categorical = adata.obs[groupby]

    obs_tidy.set_index(categorical, groupby, inplace=True)
    if gene_symbols is not None:
        # translate the column names to the symbol names
        obs_tidy.rename(
            columns=dict([(var_names[x], symbols[x]) for x in range(len(var_names))]),
            inplace=True,
        )
    categories = obs_tidy.index.categories

    return categories, obs_tidy

def stacked_violin_t(
    adata: AnnData,
    var_names: Union[_VarNames, Mapping[str, _VarNames]],
    groupby: Optional[str] = None,
    log: bool = False,
    use_raw: Optional[bool] = None,
    num_categories: int = 7,
    figsize: Optional[Tuple[float, float]] = None,
    dendrogram: Union[bool, str] = False,
    gene_symbols: Optional[str] = None,
    var_group_positions: Optional[Sequence[Tuple[int, int]]] = None,
    var_group_labels: Optional[Sequence[str]] = None,
    standard_scale: Optional[Literal['var', 'obs']] = None,
    var_group_rotation: Optional[float] = None,
    layer: Optional[str] = None,
    stripplot: bool = False,
    jitter: Union[float, bool] = False,
    size: int = 1,
    scale: Literal['area', 'count', 'width'] = 'width',
    order: Optional[Sequence[str]] = None,
    show: Optional[bool] = None,
    save: Union[bool, str, None] = None,
    row_palette: str = 'muted',
    **kwds,
):
    """\
    Stacked violin plots.
    Makes a compact image composed of individual violin plots
    (from :func:`~seaborn.violinplot`) stacked on top of each other.
    Useful to visualize gene expression per cluster.
    Wraps :func:`seaborn.violinplot` for :class:`~anndata.AnnData`.
    Parameters
    ----------
    {common_plot_args}
    stripplot
        Add a stripplot on top of the violin plot.
        See :func:`~seaborn.stripplot`.
    jitter
        Add jitter to the stripplot (only when stripplot is True)
        See :func:`~seaborn.stripplot`.
    size
        Size of the jitter points.
    order
        Order in which to show the categories. Note: if `dendrogram=True`
        the categories order will be given by the dendrogram and `order`
        will be ignored.
    scale
        The method used to scale the width of each violin.
        If 'width' (the default), each violin will have the same width.
        If 'area', each violin will have the same area.
        If 'count', a violin’s width corresponds to the number of observations.
    row_palette
        The row palette determines the colors to use for the stacked violins.
        The value should be a valid seaborn or matplotlib palette name
        (see :func:`~seaborn.color_palette`).
        Alternatively, a single color name or hex value can be passed,
        e.g. `'red'` or `'#cc33ff'`.
    standard_scale
        Whether or not to standardize a dimension between 0 and 1,
        meaning for each variable or observation,
        subtract the minimum and divide each by its maximum.
    swap_axes
         By default, the x axis contains `var_names` (e.g. genes) and the y axis the `groupby` categories.
         By setting `swap_axes` then x are the `groupby` categories and y the `var_names`. When swapping
         axes var_group_positions are no longer used
    {show_save_ax}
    **kwds
        Are passed to :func:`~seaborn.violinplot`.
    Returns
    -------
    List of :class:`~matplotlib.axes.Axes`
    Examples
    -------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
    >>> sc.pl.stacked_violin(adata, markers, groupby='bulk_labels', dendrogram=True)
    Using var_names as dict:
    >>> markers = {{'T-cell': 'CD3D', 'B-cell': 'CD79A', 'myeloid': 'CST3'}}
    >>> sc.pl.stacked_violin(adata, markers, groupby='bulk_labels', dendrogram=True)
    See also
    --------
    rank_genes_groups_stacked_violin: to plot marker genes identified using the :func:`~scanpy.tl.rank_genes_groups` function.
    """
    import seaborn as sns  # Slow import, only import if called

    if use_raw is None and adata.raw is not None:
        use_raw = True
    var_names, var_group_labels, var_group_positions = check_var_names_type(
        var_names, var_group_labels, var_group_positions
    )
    has_var_groups = (
        True
        if var_group_positions is not None and len(var_group_positions) > 0
        else False
    )
    categories, obs_tidy = prepare_dataframe(
        adata,
        var_names,
        groupby,
        use_raw,
        log,
        num_categories,
        gene_symbols=gene_symbols,
        layer=layer,
    )

    if standard_scale == 'obs':
        obs_tidy = obs_tidy.sub(obs_tidy.min(1), axis=0)
        obs_tidy = obs_tidy.div(obs_tidy.max(1), axis=0).fillna(0)
    elif standard_scale == 'var':
        obs_tidy -= obs_tidy.min(0)
        obs_tidy = (obs_tidy / obs_tidy.max(0)).fillna(0)
    elif standard_scale is None:
        pass
    else:
        logg.warning('Unknown type for standard_scale, ignored')

    if 'color' in kwds:
        row_palette = kwds['color']
        # remove color from kwds in case is set to avoid an error caused by
        # double parameters
        del kwds['color']
    if 'linewidth' not in kwds:
        # for the tiny violin plots used, is best
        # to use a thin lindwidth.
        kwds['linewidth'] = 0.5
    
    #if swap_axes:
    # plot image in which x = group by and y = var_names
    dendro_width = 3 if dendrogram else 0
    vargroups_height = 0.45 if has_var_groups else 0
    if figsize is None:
        width = len(var_names) * 0.3 + dendro_width
        height = len(categories) * 0.4 + vargroups_height
    else:
        width, height = figsize

    fig = pl.figure(figsize=(width, height))

    # define a layout of nrows = 1 x var_names columns
    # if plot dendrogram a col is added
    # each col is one violin plot.
    num_cols = len(var_names) + 1  # +1 to account for dendrogram
    width_ratios = [dendro_width] + ([1] * len(var_names))

    axs = gridspec.GridSpec(
        # 20200116 by yuanzan
        nrows=2,
        ncols=num_cols,
        height_ratios=[width - vargroups_height, vargroups_height],
        wspace=0,
        width_ratios=width_ratios,
    )

    axs_list = []
    if dendrogram:
        dendro_ax = fig.add_subplot(axs[0])
        _plot_dendrogram(
            dendro_ax, adata, groupby, orientation='top', dendrogram_key=dendrogram,
        )
        axs_list.append(dendro_ax)
    first_ax = None
    if is_color_like(row_palette):
        row_colors = [row_palette] * len(var_names)
    else:
        row_colors = sns.color_palette(row_palette, n_colors=len(var_names))
    for idx, y in enumerate(obs_tidy.columns):
        ax_idx = idx + 1  # +1 to account that idx 0 is the dendrogram
        if first_ax is None:
            ax = fig.add_subplot(axs[0, ax_idx])
            first_ax = ax
        else:
            ax = fig.add_subplot(axs[0, ax_idx])
        axs_list.append(ax)
        ax = sns.violinplot(
            x=y,
            y=obs_tidy.index,
            data=obs_tidy,
            inner=None,
            orient='horizion',
            scale=scale,
            ax=ax,
            #color=row_colors,
            **kwds,
        )

        ax.set_ylabel("")

        ax.tick_params(bottom=False, top=False, left=False, right=False)

        if idx == (len(obs_tidy.columns) -1): 
            ax.yaxis.tick_right()
        else:
            ax.set_yticklabels("")

        if stripplot:
            ax = sns.stripplot(
                # 20200116 by yuanzan
                x=obs_tidy.index,
                y=y,
                data=obs_tidy,
                jitter=jitter,
                color='black',
                size=size,
                ax=ax,
            )

        ax.set_xlabel(
                var_names[idx],
                rotation=45,
                fontsize='small',
                labelpad=8,
                ha='center',
                va='top',
        )
        ax.grid(False)
        ax.tick_params(
            axis='x',
            #right=True,
            top=False,
            #labelright=True,
            left=False,
            bottom=False,
            labelleft=False,
            labeltop=False,
            labelbottom=False,
            labelrotation=0,
            labelsize='x-small',
        )
    pl.subplots_adjust(wspace=0, hspace=0)

    _utils.savefig_or_show('stacked_violin', show=show, save=save)

    return axs_list

