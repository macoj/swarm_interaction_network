__author__ = 'marcos'
import matplotlib
import numpy
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import cycle
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram
from fastcluster import *


class Plotter:
    def __init__(self):
        pass

    @staticmethod
    def plot_curve(pd_data, title=None, x_label=None, y_label=None, output_filename=None, legends=None):
        if pd_data is not None:
            font = {'family': 'normal',
                    'weight': 'normal',
                    'size': 8}
            matplotlib.rc('font', **font)
            fig = plt.figure()
            #ax3 = fig.add_subplot(plot_gridspec[4, 2])
            if type(pd_data) is not list:
                plt.plot(pd_data['x'], pd_data['y'], linestyle='-', marker='.')
                plt.xlim(min(pd_data['x']), max(pd_data['x']))
            else:
                lines = ["-", "--", "-.", ":"]  # matplotlib.markers.MarkerStyle.markers.keys() #
                linecycler = cycle(lines)
                if legends:
                    legendcycler = cycle(legends)
                for pd_data_i in pd_data:
                    legend_title = None
                    if legends:
                        legend_title = cycle(legendcycler)
                    plt.plot(pd_data_i['x'], pd_data_i['y'], linestyle=next(linecycler),
                             marker='.', label=next(legend_title))
                    #plt.legend(loc=2)
                    plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
                plt.xlim(min(pd_data_i['x']), max(pd_data_i['x']))

            if x_label:
                plt.xlabel(x_label)
            if y_label:
                plt.ylabel(y_label)
            if title:
                plt.suptitle(title)
            if output_filename:
                plt.savefig(output_filename)
                #plt.clf()
                plt.close()
            else:
                plt.show()

    @staticmethod
    def create_heatmap(matrix,
                       main_title=None,
                       output_filename=None,
                       pd_data_1=None,
                       pd_data_2=None,
                       data_hist_1=None,
                       data_hist_2=None,
                       graphs=None,
                       ordered=False,
                       fitness_slice=500):
        matrixdf = pd.DataFrame(matrix)
        font = {'family': 'normal',
                'weight': 'normal',
                'size': 8}
        last_data_1 = 0.0
        matplotlib.rc('font', **font)
        # look at raw data
        #axi = plt.imshow(matrixdf,interpolation='nearest')
        #ax = axi.get_axes()

        #plt.clean_axis(ax)

        # row clusters
        if ordered:
            row_pairwise_dists = squareform(pdist(matrixdf))
            row_clusters = linkage(row_pairwise_dists, method='complete')
            row_dendogram = dendrogram(row_clusters, no_plot=True, count_sort='ascending')

        # calculate pairwise distances for columns
        if ordered:
            col_pairwise_dists = squareform(pdist(matrixdf.T))
            col_clusters = linkage(col_pairwise_dists, method='complete')
            col_dendogram = dendrogram(col_clusters, no_plot=True, count_sort='ascending')

        # plot the results
        fig = plt.figure(figsize=(12.5, 10))
        #plot_gridspec = gridspec.GridSpec(3,2, wspace=0.05,
        #  hspace=0.05, width_ratios=[0.25,1],height_ratios=[0.25,1,0.25])
        plot_gridspec = gridspec.GridSpec(5, 5, width_ratios=[0.15, 0.15, 0.2, 0.2, 0.2])

        ### col dendrogram ####
        #col_denAX = fig.add_subplot(plot_gridspec[0,1])
        if pd_data_1 is not None:
            title = ''
            if type(pd_data_1) == tuple:  # not so pythonic
                title = pd_data_1[0]
                pd_data_1 = pd_data_1[1]
            last_data_1 = pd_data_1['y'][len(pd_data_1)-1]

            #ax3 = fig.add_subplot(plot_gridspec[0,1])
            ax1 = plt.subplot(plot_gridspec[0, 2:])
            slice_base = max(0, len(pd_data_1) - fitness_slice)
            plt.plot(pd_data_1['x'], pd_data_1['y'], linestyle='-')
            plt.xlim(slice_base, len(pd_data_1))
            plt.title(title)
    #     else:
    #         col_denAX = fig.add_subplot(plot_gridspec[0,1])
        #create an empty graph

        ### row dendrogram ###
        ## t ODO: fix that please:
        if ordered:
            pass
            #row_denAX = fig.add_subplot(plot_gridspec[1,0])
            #row_denD = dendrogram(row_clusters, orientation='right', count_sort='ascending')
            #row_denAX.get_xaxis().set_ticks([]) # removes ticks

            #slice_base = max(0, max(pd_data_1['x']) - fitness_slice)
            #plt.plot(pd_data_1['x'], pd_data_1['y'], linestyle='-')
            #plt.xlim(slice_base, len(pd_data_1))

        if graphs is not None:
            gs_index = 0
            for title_graph in graphs:
                title, graph, graph_histogram = title_graph
                ax3 = plt.subplot(plot_gridspec[gs_index, 0])
                graph = graph.to_undirected()
                # we don't care about the weight because we already are filtering here
                nx.draw(graph, node_size=2, width=0.4, with_labels=False,
                        pos=nx.spring_layout(graph, weight=None))
                plt.title(title)
                ax3 = plt.subplot(plot_gridspec[gs_index, 1])
                # let's add the histogram, but remove all 1 values
                graph_histogram_without_one = []
                for v in graph_histogram:
                    if v != 1:
                        graph_histogram_without_one.append(v)
                print str(graph_histogram_without_one)
#                 print str(graph_histogram)
                if not graph_histogram_without_one:
                    continue
                binwidth = 1
                min_bin = numpy.min(graph_histogram_without_one)
                max_bin = numpy.max(graph_histogram_without_one)
                bins = range(min_bin, max_bin+binwidth, binwidth)
                ax3.hist(graph_histogram_without_one, bins=bins, facecolor='red', alpha=0.45)
                plt.xticks(numpy.unique(graph_histogram_without_one))
                plt.tick_params(axis='both', which='major', labelsize=5)
                plt.tick_params(axis='both', which='minor', labelsize=5)
#                plt.xticks(range(numpy.min(graph_histogram_without_one),
#                           numpy.max(graph_histogram_without_one),
#                           (numpy.min(graph_histogram_without_one) +  numpy.max(graph_histogram_without_one))/5))
                #plt.xlim(1, numpy.max(graph_histogram))
                if gs_index == 0:
                    plt.title("Components size\nhistogram")
                gs_index += 1

        ### heatmap ###
        heatmap_subplot = fig.add_subplot(plot_gridspec[1:4, 2:])

        if ordered:
            pass
            axi = heatmap_subplot.imshow(matrixdf.ix[row_dendogram['leaves'], col_dendogram['leaves']],
                                         interpolation='nearest', aspect='auto', origin='lower')
        else:
            axi = heatmap_subplot.imshow(matrixdf, interpolation='nearest', aspect='auto', origin='lower')
        # removes ticks
        heatmap_subplot.get_xaxis().set_ticks([])
        heatmap_subplot.get_yaxis().set_ticks([])
        axcolor = fig.add_axes([0.91, 0.27, 0.02, 0.45])
        plt.colorbar(axi, cax=axcolor)
        #fig.tight_layout()

        if pd_data_2 is not None:
            title = ''
            if type(pd_data_2) == tuple:  # not so pythonic
                title = pd_data_2[0]
                pd_data_2 = pd_data_2[1]
            ax3 = fig.add_subplot(plot_gridspec[4, 2])
            plt.plot(pd_data_2['x'], pd_data_2['y'], linestyle='-', marker='.')
            plt.xlim(min(pd_data_2['x']), max(pd_data_2['x']))
            #plt.ylim(0, 1.1)
            plt.title(title)

        if data_hist_1 is not None:
            title = ''
            if type(data_hist_1) == tuple:  # not so pythonic
                title = data_hist_1[0]
                data_hist_1 = data_hist_1[1]
            #binwidth = 1
            ax3 = fig.add_subplot(plot_gridspec[4, 3])
            #min_bin = numpy.min(data_hist_1)
            #max_bin = numpy.max(data_hist_1)
            #bins = range(min_bin,max_bin+binwidth,binwidth)
            ax3.hist(data_hist_1,  facecolor='blue', alpha=0.45)
            #plt.xticks(numpy.unique(data_hist_1))
            plt.tick_params(axis='both', which='major', labelsize=5)
            plt.tick_params(axis='both', which='minor', labelsize=5)
            plt.title(title)

        if data_hist_2 is not None:
            title = ''
            if type(data_hist_2) == tuple:  # not so pythonic
                title = data_hist_2[0]
                data_hist_2 = data_hist_2[1]
            if data_hist_2:
                ax3 = fig.add_subplot(plot_gridspec[4, 4])
                #bins = range(min_bin,max_bin+binwidth,binwidth)
                ax3.hist(data_hist_2,  facecolor='blue', alpha=0.45)
                #plt.xticks(numpy.unique(data_hist_1))
                plt.tick_params(axis='both', which='major', labelsize=5)
                plt.tick_params(axis='both', which='minor', labelsize=5)
                plt.title(title)

        if main_title:
            if pd_data_1 is not None:
                main_title = main_title + '\n(' + str(last_data_1).strip() + ')'
            plt.suptitle(main_title)
        if output_filename:
            plt.savefig(output_filename)
            #plt.clf()
            plt.close()
        else:
            plt.show()