import sys
import numpy as np
import matplotlib
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
if (len(sys.argv) != 2):
    print ("<all_weights>")
    exit(1)

fin = open(sys.argv[1]);
Pattern_good = "good "
lst_good = []
Pattern_bad = "bad "
lst_bad = []
good_size = 1
bad_size = 1
for line in fin:
    if "good pi" in line:
        lst = line.split("good pi")[1].strip().split(" ")
        good_size = float(lst[0])
        bad_size = float(lst[-1])
    elif not "good pi" in line and (Pattern_good in line or Pattern_bad in line):
        Pattern = Pattern_good
        lst = lst_good
        if Pattern_bad in line:
            Pattern = Pattern_bad
            lst = lst_bad
        text = line.split(Pattern)[1].strip()
        params = text.split(" ")
        if (float(params[0]) >= 0 and float(params[-1]) >=0 and
            float(params[-1]) != float('nan')):
            lst.append(float(params[0]))

def commulate( com, lst, begin, t):
    i = begin
    while i < len(lst) and lst[i] < t:
        com += 1
        i += 1
    return (i, com)
matplotlib.rcParams.update({'font.size': 18})
def plot_lines(ax, name, lst_g, lst_b, min_x, max_x, nbins):
    lst1 = sorted([x for x in lst_g if x > 0.0], key = lambda x:x)
    lst2 = sorted([x for x in lst_b if x > 0.0], key = lambda x:x)
    params = []
    i = min_x
    while i < max_x:
        params.append(i)
        i += float((max_x - min_x))/nbins
    good = []
    bad = []
    begin = 0
    com = 0
    begin2 = 0
    com2 = 0
    for param in params:
        (begin, com) = commulate( com, lst1, begin, param)
        if len(lst1) == 0:
            good.append(0)
        else:
            good.append(float(com) / float(len(lst1)))
        (begin2, com2) = commulate(com2, lst2, begin2, param)
        if len(lst2) == 0:
            bad.append(0)
        else:
            bad.append(1- float(com2)/ float(len(lst2)))
    ax.plot(params, good, params, bad)
    plt.xlabel(u"\u03A8", fontsize=24)
    plt.ylabel('FP/FN rate', fontsize = 24)
   # a.xaxis.label.set_fontsize(40)
    #a.ylabel.label.set_fontsize(40)
    #ax.set_title(name)
"""max_range = 50
count_bins = 200
min_x_norm = 1.0
h = np.histogram([x[1] for x in lst_good if x[1] > min_x_norm], bins = count_bins,  range = (min_x_norm,max_range))
widths = np.diff(h[1])
good_l = len([x[1] for x in lst_good if x[1] > min_x_norm])
plt.bar(h[1][:-1], h[0]/float(good_l), widths, alpha = 0.5, color = "red")
hist, bins = np.histogram([x[1] for x in lst_bad if x[1] > min_x_norm], bins = count_bins, range = (min_x_norm,max_range))
widths = np.diff(bins)
bad_l = len([x[1] for x in lst_bad if x[1] > min_x_norm])
plt.bar(bins[:-1], hist/float(bad_l), widths, alpha = 0.5, color = "blue")
plt.show()
"""
f, axarr = plt.subplots(1, 1)
plot_lines(axarr, "plot", lst_good, lst_bad, 0.00, 0.8, 5000)
#plot_lines(axarr[0, 1], "w", 1, lst_good, lst_bad, 0, 10, 200)
#plot_lines(axarr[0, 2], "(w / ideal_pi) / min(pi1, pi2)", 3, lst_good, lst_bad,
#0.000001, 0.05, 500)
#plot_lines(axarr[1, 0], "(w / ideal_pi) / min(cov1, cov2)", 5, lst_good, lst_bad, 0.0, 0.001, 500)
#plot_lines(axarr[1, 1], "(w / ideal_pi) / min(pi_norm1, pi_norm2)", 2, lst_good,
#lst_bad, 0.0, 0.05, 500)
#plot_lines(axarr[1, 2], "test1 / min(pi_norm1_aver, pi_norm2_aver)", 6, lst_good, lst_bad, 0.0,  3, 500)
plt.show()
