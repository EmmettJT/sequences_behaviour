import scipy.io
import os, importlib
import matplotlib.pyplot as plt
import statistics
import scipy.stats
import numpy as np
import pandas as pd
from ast import literal_eval
import pickle
import matplotlib.patches as mpatches

def conactinate_nth_items(startlist):
    concatinated_column_vectors = []
    for c in range(len(max(startlist, key=len))):
        column = []
        for t in range(len(startlist)):
            if c <= len(startlist[t])-1:
                column = column + [startlist[t][c]]
        concatinated_column_vectors.append(column)
    return concatinated_column_vectors

def convolve_movmean(y,N):
    y_padded = np.pad(y, (N//2, N-1-N//2), mode='edge')
    y_smooth = np.convolve(y_padded, np.ones((N,))/N, mode='valid') 
    return y_smooth

def plot_across_animals(data,xlabel,ylabel,Animal_ID,windowsize):

    fig,ax = plt.subplots(1, 1, figsize=(20, 10))
    ax.set_ylim([0, 1])
    for animals in range(len(Animal_ID)):
        ax.plot(data[animals][0:int(len(data[animals])*0.8)], color = 'grey',alpha = 0.5) ## only plot 90% of the data to avoid end of last session dip dominating the moving average. 

    cc_vecs = conactinate_nth_items(data)   
    medianPPerfs = [] 
    for item in cc_vecs:
        medianPPerfs = medianPPerfs + [np.median(item)]
    covlved_medianPPerfs = convolve_movmean(medianPPerfs,windowsize)

    ax.plot(covlved_medianPPerfs[0:int(len(covlved_medianPPerfs)*0.8)], color = 'firebrick')  

    patch1 = mpatches.Patch(color='grey', label=('n = ' + str(len(Animal_ID))))
    plt.legend(handles=[patch1])

    ax.set_xlabel(xlabel,fontsize = 20)
    ax.set_ylabel(ylabel,fontsize = 20)
    
    ax.set_xlim(0,8000)

def SaveFig(file_name,figure_dir):
    if not os.path.isdir(figure_dir):
        os.makedirs(figure_dir)
    plt.savefig(figure_dir + file_name, bbox_inches='tight')
    plt.close()
    
def mean_learning_curve(AA_PerfectScores,Animal_ID):

    TrialbyTrial_Pscores = conactinate_nth_items(list(np.array(AA_PerfectScores)))
    MeanLearningCurve = []
    for i,item in enumerate(TrialbyTrial_Pscores):
        MeanLearningCurve = MeanLearningCurve + [np.mean(item)]

    return MeanLearningCurve