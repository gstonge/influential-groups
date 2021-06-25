#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module to create joyplots for marginal distribution

Author: Guillaume St-Onge <guillaume.st-onge.4@ulaval.ca>
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from scipy.interpolate import interp1d
from scipy import interpolate

#figure params
fontsize=8
plt.style.use('seaborn-paper')
plt.rc('text', usetex=True)
plt.rc('font',family='serif',serif='Computer Modern')
plt.rc('xtick', labelsize=fontsize)
plt.rc('ytick', labelsize=fontsize)
plt.rc('axes', labelsize=fontsize)
color_list = [(0.816,0.820,0.902),(0.455,0.663,0.812),(0.016,0.353,0.553)]
color_map = LinearSegmentedColormap.from_list(
        'cmap', color_list)#cm.get_cmap('bone')
lw = 1.5

def smooth_curve(xvalues,yvalues,xvec,std):
    m1 = np.outer(xvec,np.ones(len(xvalues)))
    m2 = np.outer(np.ones(len(xvec)),xvalues)
    m3 = np.exp(-(m1 - m2)**2/(2*std))
    yvec = np.dot(m3,yvalues)
    return yvec

def joyplot(axes,marginal_dict, xpoints=1000, xmin=None, xmax=None, xticks=None,
            smooth=False, std=2, save=None, rescale=True, color_map=color_map,
            labelpos=(-0.28,0.02)):


    #get global xrange
    if (xmin is None) or (xmax is None):
        for marginal in marginal_dict.values():
            xvalues = list(marginal.keys())
            min_value = min(xvalues)
            max_value = max(xvalues)
            if xmin is None or min_value < xmin:
                xmin = min_value
            if xmax is None or max_value > xmax:
                xmax = max_value
    xvec = np.linspace(xmin,xmax,xpoints)

    #set axis format
    for ax in axes:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_linewidth(lw)
        ax.patch.set_facecolor('none')
        ax.set_ylim([0.,1.])
        ax.set_yticks([])
        ax.tick_params(axis=u'both', which=u'both',length=0)
        ax.set_xlim([0,xmax])
        if xticks is not None:
            ax.set_xticks(xticks)

    #order label according to marginal mean
    label_list = []
    mean_list = []
    for label in marginal_dict:
        marginal = marginal_dict[label]
        xvalues = list(marginal.keys())
        yvalues = [marginal[x]*1 for x in xvalues]
        label_list.append(label)
        mean_list.append(sum([x*y/sum(yvalues) for x,y in
                              zip(xvalues,yvalues)]))
    label_list = [label for _,label in sorted(zip(mean_list,
                                                  label_list))]

    #plot the curves
    for count,label in enumerate(reversed(label_list)):
        marginal = marginal_dict[label]
        xvalues = list(marginal.keys())
        yvalues = np.array([marginal[x]*1 for x in xvalues])

        #get mean value to specify color
        mean = sum([x*y/sum(yvalues) for x,y in zip(xvalues,yvalues)])
        # color=color_map(0.98*(mean-xmin)/(xmax-xmin) + 0.01)
        color=color_map[count]
        axes[count].spines['bottom'].set_color(color)

        if smooth==True:
            yvec = smooth_curve(xvalues,yvalues,xvec,std)
            xvalues = xvec
            yvalues = yvec.T

        if rescale:
            yvalues /= np.max(yvalues)

        #spline cubic interpolation for the look
        # f = interp1d(xvalues, 0.95*yvalues, kind='cubic')
        tck = interpolate.splrep(xvalues, 0.9*yvalues, s=0)
        xnew = np.linspace(min(xvalues),max(xvalues),1000)
        ynew = interpolate.splev(xnew,tck,der=0)

        #plot
        axes[count].plot(xnew,ynew, color='white',lw=lw)
        axes[count].fill_between(xnew,ynew,color=color)
        axes[count].text(labelpos[0],labelpos[1],"{}".format(label),
                         transform=axes[count].transAxes,
                         fontsize=fontsize,color='black')

    return axes

if __name__ == '__main__':
    marginal_dict = dict()
    for i in range(10):
        marginal_dict[i] = dict()
        for t in np.linspace(-3+i,3+i,4):
            marginal_dict[i][t] = np.exp(-(t-i)**2)

    width = 7.057/2
    nb_dist = len(marginal_dict)
    height = nb_dist*width/8
    fig, axes = plt.subplots(nb_dist,1, figsize=(width, height), sharex=True)
    plt.subplots_adjust(left=0.1, bottom=0.12, right=0.96, top=0.99,
           wspace=0.25, hspace=-0.8)

    joyplot(axes,marginal_dict, xmin=-5, xmax=15, smooth=True, std=2)

    plt.show()
