import argparse
import csv
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist

def get_nplots(nsamples, spp):
    """Return number of plots given number of samples and samples per plot (spp)"""
    nplots = 0
    while nsamples > 0:
        nsamples -= spp
        nplots += 1
    return nplots

def get_linkage(cns):
    df_dist = pdist(cns, "euclidean")
    df_link = sch.linkage(df_dist, method='average')
    return df_link

def get_cns_fclust(df_link, fclust_threshold):
    return sch.fcluster(df_link, fclust_threshold, "distance")

def reorder_df_samples(df, cns, df_link):
    df_leaves = sch.dendrogram(df_link, labels = cns.index, no_plot=True, distance_sort="descending")['leaves']
    sample_order = [cns.index[i] for i in df_leaves]
    df_new = df[["chr", "start", "end", "name"] + sample_order]
    return df_new

def plot_heatmap(df, pop_info, output_filename, color_column, hclust, annotate_column, label_heatmap, df_link, first_index = "pop", second_index = "super_pop", sample_names = True, sample_range = [0, 0],  xspace = 0.05, yspace = 0.07, xmin = 0.04, include_coords = False, fclust_threshold = None):
    """
    Plot the given DataFrame as a heatmap.
    """
    if color_column is not None:
        color_indivs = True

    col_labels = []
    color_order = ["red", "y", "g", "b", "purple", "k", "orange", "c", "firebrick", "gray", "pink", "aqua", "darkgray", "indigo", "gold", "lime", "magenta", "maroon", "olive"]
    unique_pops = []

    if color_column is None:
        pop_name, super_pop_name = "", ""
        for col in df.columns[4:]:
            out_label = ""
            if pop_info.loc[col][first_index] != pop_name:
                pop_name = pop_info.loc[col][first_index]
                out_label = pop_name
            if second_index is not None:
                if pop_info.loc[col][second_index] != super_pop_name:
                    super_pop_name = pop_info.loc[col][second_index]
                    out_label += "." + super_pop_name
            col_labels.append(out_label)
    else:
        col_colors = {}
        color_dict = {group: color_order[i] for i, group in enumerate(pop_info[color_column].unique())}
        for col in df.columns[4:]:
            pop_name = pop_info.loc[col][color_column]
            if annotate_column is not None:
                col_label = col + "." + pop_info.loc[col][annotate_column]
            else:
                col_label = col
            col_labels.append(col_label)
            if pop_name not in unique_pops:
                unique_pops.append(pop_name)
            col_colors[col_label] = color_dict[pop_name]

    if not sample_names:
        col_labels = ["|" for x in col_labels]

    if include_coords:
        region_labels = ["{region_name}\n{chr}:{coords:.1f}".format(region_name = df.index[x], chr = df.iloc[x]['chr'], coords = df.iloc[x]['start']/1000000.) for x in range(df.shape[0])]
    else:
        region_labels = df.index

    #colors = mpl.colors.ListedColormap(['w', 'gray', 'k', (0, 0, 0.5), (0, 0, 1), 'c', 'g', 'y', 'orange', 'red', 'firebrick'], name='cp_colormap')
    color_set = list(plt.cm.Greys(np.linspace(0,1,3))) + list(plt.cm.jet(np.linspace(0,1,8)))
    colors = mpl.colors.ListedColormap(color_set, name='cp_colormap')
    bounds = [x - 0.5 for x in range(12)]
    norm = mpl.colors.BoundaryNorm(bounds, colors.N)

    ### Make dendrogram based on example: http://nbviewer.ipython.org/github/OxanaSachenkova/hclust-python/blob/master/hclust.ipynb ###

    fig = plt.figure(figsize=(24,15))

    # Set size for subpanels
    xmin, ymin = xmin, 0.05
    hmap_w, hmap_h = 0.8, 0.6
    dendro_h = 0.25
    cbar_w = 0.02
    legend_w = 0.1
    xspace, yspace = xspace, yspace

    # [xmin, ymin, width, height]
    heatmap_dims = [xmin, ymin, hmap_w, hmap_h]
    sample_dendro_dims = [xmin, ymin + hmap_h + yspace, hmap_w, dendro_h]
    colorbar_dims = [xmin + hmap_w + xspace, ymin, cbar_w, hmap_h]
    sample_legend_dims = [xmin + hmap_w + xspace, ymin + hmap_h + yspace, cbar_w, dendro_h]

    # Plot legend
    legend_axis = fig.add_axes(sample_legend_dims)
    legend_axis.set_axis_off()
    legend_entries = []
    for pop in unique_pops:
        pop_entry = mpl.patches.Patch(color=color_dict[pop], label=pop)
        legend_entries.append(pop_entry)

    plt.legend(handles=legend_entries, bbox_to_anchor=sample_legend_dims, borderaxespad=0.)

    
    ### Unfinished code for clustering by region
    cluster_regions = False

    if cluster_regions:
        cns = df[df.columns[4:]].T
        df_dist = pdist(cns)
        Y = linkage(df_dist, method='average')
        region_dendro = fig.add_axes([0.05,0.1,0.2,0.6])
        Z1 = sch.dendrogram(Y, orientation='right',labels=cns.index) # adding/removing the axes
        ax1.set_yticks([])
        ax1.set_xticks([])
    ###


    # Compute and plot sample dendrogram.
    if hclust:
        sample_dendro = fig.add_axes(sample_dendro_dims)
        Z2 = sch.dendrogram(df_link, labels=df.columns[4:], distance_sort="descending")
        sample_dendro.set_yticks([])
        if not label_heatmap:
            sample_dendro.set_xticklabels(col_labels, rotation=90, fontsize=4)
        else:
            sample_dendro.set_xticklabels([" " for label in col_labels], fontsize=4)

        # Add grey box to indicate included samples if plotting a subset of samples
        if sample_range[0] != 0 or sample_range[1] != len(col_labels):
            xlocs = sample_dendro.xaxis.get_majorticklocs()
            xmin = xlocs[sample_range[0]]
            xmax = xlocs[sample_range[1] - 1]
            ymax = sample_dendro.get_ylim()[1]
            sample_dendro.add_patch(mpl.patches.Rectangle((xmin, 0), xmax - xmin, ymax, facecolor="grey", edgecolor="none"))

        if fclust_threshold is not None:
            min_x, max_x = sample_dendro.get_xlim()
            ymax = sample_dendro.get_ylim()[1]
            max_y = max(map(lambda x: x[1], Z2['dcoord']))
            plt.axhline(y = fclust_threshold, xmin = min_x, xmax = max_x, color="k")

    #Compute and plot the heatmap
    axmatrix = fig.add_axes(heatmap_dims)
    heatmap = axmatrix.imshow(df[df.columns[sample_range[0] + 4:sample_range[1] + 4]], aspect='auto', origin='lower', norm=norm, cmap=colors, interpolation='nearest')
    axmatrix.set_yticks(range(len(df.index)))
    axmatrix.set_yticklabels(region_labels, fontsize=10)
    axmatrix.set_xticks(range(sample_range[1] - sample_range[0]))
    axmatrix.xaxis.tick_top()
    if sample_names:
        if label_heatmap:
            axmatrix.set_xticklabels(col_labels[sample_range[0]:sample_range[1]], rotation=90, fontsize=4)
        else:
            axmatrix.set_xticklabels([])

    map(lambda x: x.set_visible(False), axmatrix.xaxis.get_majorticklines())
    map(lambda x: x.set_visible(False), axmatrix.yaxis.get_majorticklines())

    # Plot colorbar.
    axcolor = fig.add_axes(colorbar_dims)
    cbar = plt.colorbar(heatmap, cax=axcolor, ticks=range(11))
    cbar.ax.set_yticklabels([str(x) for x in range(10)] + ['10+'])

    ######
    if color_column is not None:
        if label_heatmap:
            map(lambda x: x.set_color(col_colors[x.get_text()]), axmatrix.get_xticklabels())
        else:
            map(lambda x: x.set_color(col_colors[x.get_text()]), sample_dendro.get_xticklabels())

    plt.savefig(output_filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("pop_file")
    parser.add_argument("output_file_prefix")
    parser.add_argument("--plot_type", default="png", choices=["pdf","png"], help = "Plot file format")
    parser.add_argument("--spp", default=None, type=int, help="Number of samples per plot (Default: all)")
    parser.add_argument("--color_column", default=None, help="Label every sample and color by specified column from pop file (e.g. super_pop)")
    parser.add_argument("--hclust", action="store_true", help="Group samples using hierarchical clustering")
    parser.add_argument("--fclust_threshold", type=float, default=None, help="Threshold for dividing samples into flat clusters (Default: %(default)s)")
    parser.add_argument("--exclude_sample_names", action="store_true", help="Use '-' instead of sample name")
    parser.add_argument("--annotate_column", default=None, help="Name of column to append to sample names")
    parser.add_argument("--label_heatmap", action="store_true", help="Label heatmap instead of dendrogram")
    parser.add_argument("--xspace", type=float, default = 0.05, help="Amount of space between heatmap and colorbar/legend (Default: %(default)s)")
    parser.add_argument("--yspace", type=float, default = 0.07, help="Amount of space between heatmap and dendrogram (Default: %(default)s)")
    parser.add_argument("--xmin", type=float, default = 0.04, help="Amount of whitespace left of heatmap and dendrogram (Default: %(default)s)")
    parser.add_argument("--include_coords", action="store_true", help="Annotate regions with genomic coordinates")
    args = parser.parse_args()

    df = pd.read_table(args.input_file, na_values="NA").dropna()
    pop_info = pd.read_table(args.pop_file)
    pop_info.sort(columns=["super_pop", "pop"], inplace=True, axis=0)

    sample_order = [sample for sample in pop_info.sample if sample in df.columns]

    if args.spp is None or args.spp > len(sample_order):
        spp = len(sample_order)
    else:
        spp = args.spp

    pop_info.index = pop_info.sample
    pop_info = pop_info.ix[sample_order, :]

    df = df[["chr", "start", "end", "name"] + sample_order]
    df["name"] = df.name.map(lambda x: "_".join(sorted(x.split(","))))
    df.sort(columns=["name"], inplace=True, axis=0)
    df.index = df.name
    df = df[sample_order].applymap(lambda x: 10 if x > 10 else x)

    sample_names = not args.exclude_sample_names

    # Plot heatmap
    nsamples = len(sample_order)
    nplots = get_nplots(nsamples, spp)

    if args.hclust:
        cns = df[df.columns[4:]].T
        df_link = get_linkage(cns)
        cns_fclust = get_cns_fclust(df_link, args.fclust_threshold)
        cns_fclust = pd.DataFrame(data={"sample": cns.index, "cluster": cns_fclust})
        cns_fclust.sort(inplace=True, columns="cluster")
        cns_fclust.to_csv(args.output_file_prefix + ".tab", index=False, sep="\t")
        df = reorder_df_samples(df, cns, df_link)
    else:
        df_link = None

    for i in range(nplots):
        start_sample = i * spp
        end_sample = min(i * spp + spp, nsamples)
        plot_name = ".".join([args.output_file_prefix, str(i), args.plot_type])
        plot_heatmap(df, pop_info, plot_name, args.color_column, args.hclust, args.annotate_column, args.label_heatmap, df_link, sample_names = sample_names, sample_range = [start_sample, end_sample], xspace = args.xspace, yspace = args.yspace, xmin = args.xmin, include_coords = args.include_coords, fclust_threshold = args.fclust_threshold)

