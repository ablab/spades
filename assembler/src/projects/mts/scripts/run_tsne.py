__author__ = 'tanunia'


import copy
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import bhtsne
import matplotlib.patches as mpatches
import matplotlib.cm as cm

def draw_points(points, names, fig):
    ax = fig.add_subplot(111)
    for i in xrange(len(points)):
        ax.annotate(names[i], xy=points[i], textcoords='data')

def points_mean(x, y):
    x_mean = sum(x)*1.0/len(x)
    y_mean = sum(y)*1.0/len(y)

    return x_mean, y_mean

def find_cluster_centers(mp):
    points = []
    names = []
    for c in mp.keys():
        names.append("Cl-" + str(c))
        x, y = points_mean(mp[c]['x'], mp[c]['y'])
        points.append([x, y])
    return points, names

def divide_by_cluster(names, clusters):
    res = {}
    for i in xrange(len(clusters)):
        c = clusters[i]
        if not c in res.keys():
            res[c] = []
        res[c].append(names[i])

    return res

def take_first_per(clusters, per = 0.1):
    res = []
    for c in clusters.keys():
        for i in xrange(max(int(len(clusters[c])*0.1), min(10, len(clusters[c])) ) ):
            res.append(clusters[c][i])

    return res



def divide_by_color(x, y, color):
    res = {}

    for i in xrange(len(color)):
        c = color[i]
        if not c in res.keys():
            res[c] = {}
            res[c]["x"] = []
            res[c]["y"] = []
        res[c]["x"].append(x[i])
        res[c]["y"].append(y[i])

    return res

def form_points(df):
    x = df["x"].tolist()
    y = df["y"].tolist()
    names = [z[8:18] for z in df[0].tolist()]
    points = zip(x, y)

    return points, names

import re
extract_num = re.compile("\d+")

def run_tsne(features_file, colors_file, output_prefix
             , filter_sample=[]
             , filter_cluster=[]
             , lst=[]
             , draw_per = 1.0
             , iter = 1000
             , perplexity = 50):
    # read data
    data_df = pd.read_table(features_file, header=None)
    cluster_colors = pd.read_table(colors_file, header=None)
    print(data_df.head())

    # make dataframe pretty
    cluster_colors = cluster_colors.rename(columns={1:'color'})
    cluster_colors["color"] = [int(extract_num.findall(str(x))[0]) for x in cluster_colors["color"].tolist()]
    print(cluster_colors.head())
    #cluster_colors = cluster_colors.rename(columns={0:0})

    # filter by samples
    if len(filter_sample) > 0:
        filter1 = []
        for x in cluster_colors[0].tolist():
            for it in filter_sample:
                st = "sample" + it + "-"
                if x.startswith(st):
                    filter1.append(x)
        cluster_colors = cluster_colors[cluster_colors[0].isin(filter1)]

    # filter by percent
    if draw_per < 1:
        clusters = divide_by_cluster(cluster_colors[0].tolist(), cluster_colors["color"].tolist())
        filter2 = take_first_per(clusters, lst)
        s = set(filter2)
        lst_new = []
        for n in lst:
            for x in cluster_colors[0].tolist():
                if x.startswith(n):
                    print x
                    lst_new.append(x)
                    if x not in s:
                        filter2.append(x)
        lst = lst_new
        cluster_colors = cluster_colors[cluster_colors[0].isin(filter2)]


    # merge data
    mapped = pd.merge(cluster_colors, data_df, on=0)

    # filter by length
    mapped["length"] = [int(x.split("_")[3]) for x in mapped[0].tolist()]
    mapped = mapped[mapped["length"] > 2000]
    print(mapped)

    # normalize like in CONCOCT
    data = mapped.as_matrix(columns=mapped.columns[2:-1])

    v = (1.0/mapped["length"]).as_matrix()[:, np.newaxis]
    data = data + v
    along_Y = np.apply_along_axis(sum, 0, data)
    data = data/along_Y[None, :]
    along_X = np.apply_along_axis(sum, 1, data)
    data = data/along_X[:, None]
    data = np.log(data)
    #print(data)

    embedding_array = bhtsne.run_bh_tsne(data, initial_dims=data.shape[1], perplexity=perplexity, max_iter=iter)
    mapped["x"] = embedding_array[:, 0]
    mapped["y"] = embedding_array[:, 1]

    # draw result of TSNE on scatter plot

    pp = PdfPages(output_prefix)


    # filter clusters to show
    fc = filter_cluster
    if len(fc) > 0:
        filtered = mapped[mapped["color"].isin(fc)]
        #mapped = filtered
    else:
        filtered = mapped

    fig = pyplot.figure()

    # draw scatter plot
    color = mapped["color"].tolist()
    mx_color = max(color)
    pyplot.scatter(mapped["x"].tolist(), mapped["y"].tolist(), c=[cm.spectral(float(i) /mx_color) for i in color])

    # make a legend for specific clusters
    # find cluster centers
    x = filtered["x"].tolist()
    y = filtered["y"].tolist()
    mp = divide_by_color(x, y, filtered["color"].tolist())
    points, names = find_cluster_centers(mp)
    patches = []
    dcolors = list(set(color))
    for c in dcolors:
        if c in fc and len(fc) < 5:
            patches.append(mpatches.Patch(color=cm.spectral(float(c)/mx_color), label='C-'+ str(c)))
    pyplot.legend(handles=patches)
    draw_points(points, names, fig)

    # mark specific points
    filtered = mapped[mapped[0].isin(lst)]
    pyplot.scatter(filtered["x"].tolist(), filtered["y"].tolist(), marker="p", edgecolors='black', c=[cm.spectral(float(i) /mx_color) for i in filtered["color"].tolist()])


    pyplot.title('Perp = '+ str(perplexity)+ ' Iter = ' + str(iter))
    pp.savefig()

    pp.close()

def get_points(file):
    if file == "":
        return []
    else:
        points = []
        fin = open(file, "r")
        for l in fin.readlines():
            points.append(l.strip())
        fin.close()
        return points

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("profile", help="profile information (depth)")
    parser.add_argument("binning", help="file with binning results")
    parser.add_argument("output", help="path to pdf-file to save graph")
    parser.add_argument("-p", "--percent", help="sets size of random subsample from profile to run TSNE",
                        type=float,
                        default=1.0)
    parser.add_argument("-i", "--iteration", help="number of TSNE iterations",
                        type=int,
                        default=1000)
    parser.add_argument("-e", "--perplexity", help="TSNE perplexity",
                        type=float,
                        default=50)
    parser.add_argument("-s", "--samples", help="run TSNE only on samples from the list",
                        nargs='+',
                        default=[])
    parser.add_argument("-c", "--clusters", help="draw only clusters from the list",
                        nargs='+',
                        default=[])
    parser.add_argument("-f", "--pointsfile", help="highlight specific points on the graph",
                        default="")

    args = parser.parse_args()
    points = get_points(args.pointsfile)
    run_tsne(args.profile, args.binning, args.output
             , args.samples
             , args.clusters
             , points
             , args.percent
             , args.iteration
             , args.perplexity)

if __name__ == "__main__":
    main()
