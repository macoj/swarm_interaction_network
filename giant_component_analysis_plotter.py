__author__ = 'marcos'
import matplotlib
import numpy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import cycle


class GiantComponentDeathPlotter:
    def __init__(self):
        pass

    @staticmethod
    def giant_component_death_curve(calculate_on, pd_datas, windows_size, xlim=None, ylim=None, figsize=None, tight_layout=None):
        # font = {'family': 'normal',
        #         'weight': 'normal',
        #         'size': 8}
        # matplotlib.rc('font', **font)
        # plot the results
        if figsize:
            fig = plt.figure(figsize=figsize)
        else:
            fig = plt.figure()
        plot_gridspec = gridspec.GridSpec(3, 3, width_ratios=[1, 0.001, 1], height_ratios=[1, 0.001, 1])
        graphs_index = 0
        ylabel = "Giant component size"
        ylabel = "Number of components"
        xlabel = "Normalized edge weight"
        markers = ['D', 's', '|', 'x', '_', '^', 'd', 'h', '+', '*', ',', 'o', '.', '1', 'p', '3', '2', '4', 'H', 'v', '8', '<', '>']
        markers = ['D', 's', 'x', '^', 'd', 'h', '+', '*', ',', 'o', '.', '1', 'p', '3', '2', '4', 'H', 'v', '8', '<', '>']
        colors = ['#e41a1c', '#377eb8', '#b2df8a', '#33a02c']
        colorblind = [(0, 0, 0), (230, 159, 0), (86, 180, 233), (0, 158, 115), (240, 228, 66), (0, 114, 178), (213, 94, 0), (204, 121, 167)]
        colorblind = ['#%02x%02x%02x' % (c[0], c[1], c[2]) for c in colorblind]
        colors = colorblind[3], colorblind[5], colorblind[6], colorblind[1]
        colorcycler = cycle(colors)
        colorcycler.next()
        x_lim = [numpy.inf, 0]
        y_lim = [numpy.inf, 0]
        for title, pd_data in pd_datas:
            for title_legend in pd_data:
                x_lim[0] = min(x_lim[0], min(pd_data[title_legend]['x']))
                x_lim[1] = max(x_lim[1], max(pd_data[title_legend]['x']))
                y_lim[0] = min(y_lim[0], min(pd_data[title_legend]['y']))
                y_lim[1] = max(y_lim[1], max(pd_data[title_legend]['y']))
        for title, pd_data in pd_datas:
            x = graphs_index / (len(pd_datas) / 2)
            y = graphs_index % (len(pd_datas) / 2)
            y = 2 if y == 1 else y
            x = 2 if x == 1 else x
            #ax3 = fig.add_subplot(plot_gridspec[x, y])
            ax = fig.add_subplot(plot_gridspec[x, y])
            ax.set_axisbelow(True)
            legends = pd_data.keys()
            legends.sort()
            # lines = ["-", "--", "-.", ":"]  # matplotlib.markers.MarkerStyle.markers.keys() #
            # markers = ["."]*len(lines) + ["v"]*len(lines) + ["^"]*len(lines)
            #
            # linecycler = cycle(lines)
            markercycler = cycle(markers)
            for title_legend in legends:
                pd_data_i = pd_data[title_legend]
                plt.grid(True)
                color = next(colorcycler)
                print 'color:' + color
                plt.plot(pd_data_i['x'], pd_data_i['y'],
                         # linestyle=next(linecycler),
                         marker=next(markercycler),
                         label=title_legend, color=color, mec=color, mew=1.4, linewidth=1.2,
                         markevery=2,
                         markersize=3)
            if not xlim:
                plt.xlim(x_lim[0], x_lim[1])
            else:
                plt.xlim(*xlim)
            if not ylim:
                plt.ylim(y_lim[0], y_lim[1])
            else:
                plt.ylim(*ylim)
            plt.title(title)
            graphs_index += 1
        plt.legend(loc=4,  numpoints=1)
        ax3 = fig.add_subplot(plot_gridspec[1, 0])
        ax3.set_yticks([])
        ax3.set_xticks([])
        ax3.set_frame_on(False)
        plt.ylabel(ylabel, fontsize=9, labelpad=20)
        ax3 = fig.add_subplot(plot_gridspec[2, 1])
        ax3.set_yticks([])
        ax3.set_xticks([])
        ax3.set_frame_on(False)

        plt.xlabel(xlabel, fontsize=9, labelpad=20)
        # plt.suptitle("Snapshot of the " + str(calculate_on) + "th iteration", fontsize=14)
        # plt.savefig('/home/marcos/PhD/research/pso_influence_graph_communities/giant' +
        #             str(calculate_on) + '-' +
        #             '_'.join(map(str, windows_size)) + '.png', bbox_inches='tight')
        # plt.close()
        if tight_layout is not None:
            if not tight_layout:
                plt.tight_layout()
            else:
                plt.tight_layout(rect=tight_layout)
        plt.show()

    @staticmethod
    def giant_component_death_curve_with_area(pd_datas, xlim=None, ylim=None, figsize=None, output_filename=None,
                                              tight_layout=None, **kargs):
        # font = {'family': 'normal',
        #         'weight': 'normal',
        #         'size': 10}
        # matplotlib.rc('font', **font)
        # plot the results
        # fig = plt.figure(figsize=(9, 6))
        # plot_gridspec = gridspec.GridSpec(2, 2, width_ratios=[1, 0.001, 1], height_ratios=[1, 0.001, 1])
        if not figsize:
            figsize = (10, 2)
        f, axs = plt.subplots(1, len(pd_datas), sharey=True, figsize=figsize)
        graphs_index = 0
        ylabel = "Number of components"
        xlabel = "Normalized edge weight"
        markers = ['D', 's', '|', 'x', '_', '^', 'd', 'h', '+', '*', ',', 'o', '.', '1', 'p', '3', '2', '4', 'H', 'v', '8', '<', '>']
        markers = ['D', 's', 'x', '^', 'd', 'h', '+', '*', ',', 'o', '.', '1', 'p', '3', '2', '4', 'H', 'v', '8', '<', '>']
        markers = ['s', 'o']
        colors = ["#e41a1c", "#377eb8", "#fdae61", "#4daf4a"]
        x_lim = [numpy.inf, 0]
        y_lim = [numpy.inf, 0]
        legends = []
        for _, pd_data in pd_datas:
            legends = pd_data.keys()
            legends.sort()
            for title_legend in pd_data:
                x_lim[0] = min(x_lim[0], min(pd_data[title_legend]['x']))
                x_lim[1] = max(x_lim[1], max(pd_data[title_legend]['x']))
                y_lim[0] = min(y_lim[0], min(pd_data[title_legend]['y']))
                y_lim[1] = max(y_lim[1], max(pd_data[title_legend]['y']))
        axs[0].set_ylabel(ylabel)
        for title_legend in legends:
            ax = axs[graphs_index]
            # lines = ["-", "--", "-.", ":"]  # matplotlib.markers.MarkerStyle.markers.keys() #
            # markers = ["."]*len(lines) + ["v"]*len(lines) + ["^"]*len(lines)
            #
            # linecycler = cycle(lines)
            markercycler = cycle(markers)
            colorcycler = cycle(colors)
            hatch = cycle(['.', '|'])
            for title, pd_data in pd_datas:
                color = next(colorcycler)
                #curve = pd_datas[title][title_legend]
                # pd_data = curve
                ax.grid(True)
                ax.plot(list(pd_data[title_legend]['x']) + [1.0], list(pd_data[title_legend]['y']) + [100],
                         # linestyle=next(linecycler),
                         marker=next(markercycler),
                         label=title, color=color, mec=color,
                         markersize=3, markevery=3, **kargs)
                ax.fill_between(list(pd_data[title_legend]['x']) + [1.0], 0.0, list(pd_data[title_legend]['y']) + [100.0], facecolor=color, alpha=0.15) #, hatch=next(hatch))
            if not xlim:
                ax.set_xlim(x_lim[0], x_lim[1])
            else:
                ax.set_xlim(*xlim)
            if not ylim:
                ax.set_ylim(y_lim[0], y_lim[1])
            else:
                ax.set_ylim(*ylim)
            ax.set_title("$t_w=%d$" % title_legend)
            graphs_index += 1
        axs[1].legend(loc=4, numpoints=1)
        # plt.suptitle("Snapshot of the " + str(calculate_on) + "th iteration", fontsize=14)
        # plt.xlabel(xlabel, fontsize=12, labelpad=20)
        f.text(0.5, 0.04, xlabel, ha='center', va='center', fontsize=9)
        plt.tight_layout()
        plt.subplots_adjust(0.09, 0.15)
        if tight_layout is not None:
            if not tight_layout:
                plt.tight_layout()
            else:
                plt.tight_layout(rect=tight_layout)
        if output_filename:
            plt.savefig(output_filename)
            #plt.clf()
            plt.close()
        else:
            plt.show()