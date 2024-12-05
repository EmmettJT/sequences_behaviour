import scipy.io
import os, importlib
import matplotlib.pyplot as plt
import statistics
import scipy.stats
import matplotlib.patches as mpatches
import seaborn as sns; sns.set()
import tkinter as tk
from matplotlib import gridspec
import matplotlib.colors
import numpy as np
import pandas as pd
import ptitprince as pt
from ast import literal_eval
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from IPython.display import HTML



def Find_Transition_times(TransitionLatency_Tfilt,TransitionTypes_Tfilt,Transitions):
    transit_times = []
    for Transits in Transitions:
        transit_times_temp = []
        for ind, transition_pair in enumerate(TransitionTypes_Tfilt):
            if transition_pair == Transits:
                transit_times_temp = transit_times_temp + [TransitionLatency_Tfilt[ind]]
                
        transit_times = transit_times + [transit_times_temp]   
    return(transit_times)


def Find_Transition_times_sns(TransitionLatency_Tfilt,TransitionTypes_Tfilt,Transitions):
    transit_times = []
    transit_type = []
    for i,Transits in enumerate(Transitions):
        transit_times_temp = []
        for ind, transition_pair in enumerate(TransitionTypes_Tfilt):
            if transition_pair == Transits:
                transit_times_temp = transit_times_temp + [TransitionLatency_Tfilt[ind]]
                
        transit_times = transit_times + transit_times_temp
        transit_type = transit_type + len(transit_times_temp)*[i+1]
    return(transit_times,transit_type)

def SaveFig(file_name,figure_dir):
    if not os.path.isdir(figure_dir):
        os.makedirs(figure_dir)
    plt.savefig(figure_dir + file_name,bbox_inches=0)
    plt.close()
    
def port_fitted_poke_times(Fitted_tfiltered_seqs,Fitted_tfiltered_times,port,max_filter):
    PortPokes = [[],[],[],[],[],[],[],[]]
    for index, seq in enumerate(Fitted_tfiltered_seqs):
        seq = literal_eval(seq) 
        if np.size(seq) > 0:
            if int(str(seq[0])[0]) == port: # if sequence starts with port
                if not int(str(seq[0])[1]) == port: #ignore self pokes 
                    current_seq_time = 0
                    for ind,item in enumerate(seq):
                        if ind > max_filter: # ignore long chains of seqs that dont return to the start port as they skew data towards being super long..
                            break
                        c_port = int(str(item)[-1])-1
                        PortPokes[c_port] = PortPokes[c_port] +[current_seq_time + literal_eval(Fitted_tfiltered_times[index])[ind]]
                        current_seq_time = current_seq_time + literal_eval(Fitted_tfiltered_times[index])[ind]
    return PortPokes


def create_plotting_df(portpokes,new_order):
    concatinated = []
    ids = []
    for index, port in enumerate(new_order):
        concatinated = concatinated + portpokes[port] 
        if index < 5:
            ids = ids +  (len(portpokes[port]) * [index])
        else:
            ids = ids +  (len(portpokes[port]) * [5])
            
    df = pd.DataFrame({'index' : ids, 
                       'Time':concatinated})
    return(df)
    
def Moving_median(Var,window):
    cumsum, moving_meds,uqs,lqs = [0], [],[],[]
    for i, x in enumerate(Var, 1):
        if i>=window:
            moving_med = np.median(Var[i-window:i])
            lst = sorted(Var[i-window:i])
            upperquart = lst[int(0.75 * window)] 
            lowerquart = lst[int(0.25 * window)]
            #can do stuff with moving_ave here
            moving_meds.append(moving_med)
            uqs.append(upperquart)
            lqs.append(lowerquart)
    return moving_meds,uqs,lqs

def convolve_movmean(y,N):
    y_padded = np.pad(y, (N//2, N-1-N//2), mode='edge')
    y_smooth = np.convolve(y_padded, np.ones((N,))/N, mode='valid') 
    return y_smooth


def mklist(n):
    for _ in range(n):
        yield []

def ten_s_filter_df(df1,filt_time):

    filtered_times = np.array(df1.loc[:,'Time'])[np.where(np.array(df1.loc[:,'Time']) < filt_time )]
    filtered_Ids = np.array(df1.loc[:,'index'])[np.where(np.array(df1.loc[:,'Time']) < filt_time )]
    
    new_df = pd.DataFrame({'index' : filtered_Ids, 
                   'Time':filtered_times})
    return(new_df)

def return_sequence_accuracy_props(TimeFiltered_seqs,Seq_references_concat):

    Correct = [21,16,63,37,72]
    Error = [23,24,25,26,27,28,12,13,14,15,17,18,61,62,64,65,67,68,31,32,34,35,36,38]
    Neutral = [11,22,33,66,41,42,43,44,45,46,47,48,51,52,53,54,55,56,57,58,71,73,74,75,76,77,78,81,82,83,84,85,86,87,88]

    Correct_scores = []
    Error_scores = []
    Neutral_scores = []
    x_refs = []
    for ind,sequence in enumerate(TimeFiltered_seqs):
        current_seq = literal_eval(sequence)
        if str(current_seq[0])[0] == '7':
            current_seq = current_seq[1::]
        if np.size(current_seq) > 0:
            Correct_scores = Correct_scores + [determine_seq_score(current_seq,Correct)]
            Error_scores = Error_scores + [determine_seq_score(current_seq,Error)]
            Neutral_scores = Neutral_scores + [determine_seq_score(current_seq,Neutral)]



    Correct_proportion = []
    Error_proportion = []
    Neutral_proportion = []
    for i in range(len(Correct_scores)):
        total = Correct_scores[i] + Error_scores[i]+ Neutral_scores[i]
        if total > 0:
            x_refs = x_refs + [Seq_references_concat[i]]
            Correct_proportion = Correct_proportion + [Correct_scores[i]/total]
            Error_proportion = Error_proportion + [Error_scores[i]/total]
            Neutral_proportion = Neutral_proportion + [Neutral_scores[i]/total]

    return Correct_proportion,Error_proportion,Neutral_proportion,x_refs

def determine_seq_score(current_seq,variable_list):
    score = 0
    for item in variable_list:
        score= score + current_seq.count(item)
    return(score)


def determineTransitionNumber(TimeFiltered_seqs):
    #Set needed reference frame:
    TransitionTypesIndex = np.array([11,12,13,14,15,16,17,18,21,22,23,24,25,26,27,28,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,51,52,53,54,55,56,57,58,61,62,63,64,65,66,67,68,71,72,73,74,75,76,77,78,81,82,83,84,85,86,87,88])
    trajects = []
    for inds, seqs in enumerate(TimeFiltered_seqs):
        seqs = literal_eval(seqs) # convert back from weird df string conversion thing
        for ind, transits in enumerate(seqs):
#             if not str(transits)[0] == str(transits)[1]:
            trajects = np.append(trajects,transits)
    transition_number = []
    for transit_types in TransitionTypesIndex:
        temp = (np.where(trajects == float(transit_types)))
        transition_number.append(len(temp[0]))
    return transition_number

############

def reversedata(port_transits,myorder):

    reordered_port_transits = []
    for i in range(1,len(port_transits)+1):
        mylist = port_transits[-i]
        mylist = [mylist[i] for i in myorder]
        newlist = []
        for item in mylist:
            newlist = newlist + [float(item)]
        reordered_port_transits = reordered_port_transits + [newlist]
    #restructure data to swap x and y axis:
    data = [[],[],[],[],[],[],[],[]]
    for ind in range(8):
        for index,item in enumerate(reordered_port_transits):
            data[ind] = data[ind] + [item[len(item)-1-ind]]
    for i in range(8):
        data[i].reverse()
    return data

def determime_heatmapdata(var,port1,port2,port3,port4):
    port1_transits = []
    for i in range(((port1*8)-8),((port1*8)-8)+8):
        port1_transits = port1_transits + [var[i]]

    port2_transits = []
    for i in range(((port2*8)-8),((port2*8)-8)+8):
        port2_transits = port2_transits + [var[i]]

    port3_transits = []
    for i in range(((port3*8)-8),((port3*8)-8)+8):
        port3_transits = port3_transits + [var[i]]

    port4_transits = []
    for i in range(((port4*8)-8),((port4*8)-8)+8):
        port4_transits = port4_transits + [var[i]]

    port_transits = [port1_transits] + [port2_transits] + [port3_transits] + [port4_transits]

    return port_transits


def plot_transition_heatmap(transition_number, port1, port2, port3, port4, port5, CurrentAnimal, file, Experiment_type, InputPathCurrent):
    # put data into correct format
    port1_transits = []
    for i in range(((port1*8)-8),((port1*8)-8)+8):
        port1_transits = port1_transits + [transition_number[i]]
    port2_transits = []
    for i in range(((port2*8)-8),((port2*8)-8)+8):
        port2_transits = port2_transits + [transition_number[i]]
    port3_transits = []
    for i in range(((port3*8)-8),((port3*8)-8)+8):
        port3_transits = port3_transits + [transition_number[i]]
    port4_transits = []
    for i in range(((port4*8)-8),((port4*8)-8)+8):
        port4_transits = port4_transits + [transition_number[i]]
    port_transits = [port1_transits] + [port2_transits] + [port3_transits] + [port4_transits]
    
    # reorder data so that it is correct for heatmap:
    a = np.array([0,1,2,3,4,5,6,7])
    a = np.delete(a, [port1-1,port2-1,port3-1,port4-1,port5-1])
    new_order = [port1-1] + [port2-1] + [port3-1] + [port4-1] + [port5-1] + list(a)
    reordered_port_transits = []
    for i in range(1,len(port_transits)+1):
        mylist = port_transits[-i]
        myorder = new_order
        mylist = [mylist[i] for i in myorder]
        newlist = []
        for item in mylist:
            newlist = newlist + [float(item)]
        reordered_port_transits = reordered_port_transits + [newlist]
    
    # restructure data to swap x and y axis:
    data = [[],[],[],[],[],[],[],[]]
    for ind in range(8):
        for index,item in enumerate(reordered_port_transits):
            data[ind] = data[ind] + [item[len(item)-1-ind]]
    for i in range(8):
        data[i].reverse()
    
    x_axis_labels = ['Port 1','Port 2','Port 3','Port 4'] # labels for x-axis
    y_axis_labels = ['Port Z','Port Y','Port X','Port 5','Port 4','Port 3','Port 2','Port 1'] # labels for y-axis
    labels =  np.array([['','','',''],
                        ['','','',''],
                        ['','','',''],
                        ['','','','T4'],
                        ['','','T3',''],
                        ['','T2','',''],
                        ['T1','','',''],
                        ['','','','']])
    mask = np.zeros_like(data)
    mask[4][3] = 1
    mask[5][2] = 1
    mask[6][1] = 1
    mask[7][0] = 1
    
    # plot:
    with sns.axes_style("white"):
        f, ax = plt.subplots(figsize=(5, 15))
        ax = sns.heatmap(data, xticklabels=x_axis_labels, yticklabels=y_axis_labels, 
                         linewidths=.5, mask=mask, square=True, cmap="YlGnBu", annot=labels, fmt='', cbar_kws=dict(use_gridspec=False, location="top"))
    
    # add separating line:
    ax.hlines([3], *ax.get_ylim())
    
    # labels:
    plt.ylabel('End Port', size=20)
    plt.xlabel('Start Port', size=20)
    fig = ax.get_figure()
    
    # save
    SaveFig((CurrentAnimal + '_' + file + '_' + Experiment_type[0][2::] + '_' + 'TransitionHeatmap.png'), InputPathCurrent + '\\SingleSessions\\' + file + '\\')



def plot_poke_proportions(transition_data, TransitionLatency_unfilt, file, Experiment_type, CurrentAnimal, InputPathCurrent, port1, P1alignedtimefiltseqs_data, new_order):
    import matplotlib.pyplot as plt

    # pull data from dataframe 
    StartPort = list(transition_data.loc[:, 'Start_Port'])   
    EndPort = list(transition_data.loc[:, 'Start_Port'])

    # determine number of pokes per port (double pokes filtered by 2s)
    filter_time = 2
    poke_count = [0] * 8
    for index, port in enumerate(StartPort):
        if port == EndPort[index]:
            if TransitionLatency_unfilt[index] < filter_time:
                poke_count[port - 1] += 1
        else:
            poke_count[port - 1] += 1

    # plot:        
    fig, ax = plt.subplots(1, 1, figsize=(20, 20))
    ax.set_ylim([0, 5])
    ax.set_xlim([0, 5])
    normalise = max(poke_count)
    circles = [plt.Circle((i % 4 + 1, i // 4 + 2), 0.5 * (poke_count[i] / normalise), fill=False, linewidth=5, color='grey') for i in range(8)]
    for circle in circles:
        ax.add_artist(circle)
    plt.axis('off')
    plt.text(0.7, 1, ('Port_Poke_proportions for Session ' + str(file) + '_' + str(Experiment_type[0])), horizontalalignment='left', size=20)
    for i in range(8):
        plt.text((i % 4 + 0.9), (i // 4 + 1.9), (str(round((poke_count[i] / sum(poke_count)) * 100, 1)) + '%'), size=20)

    # save:
    SaveFig((CurrentAnimal + '_' + file + '_' + Experiment_type[0][2::] + '_' + 'PokeProportions.png'), InputPathCurrent + '\\SingleSessions\\' + file + '\\')

    matplotlib.style.use('default')

    ## Port1 fitted poke histograms:
    # pull in data 
    P1Fitted_tfiltered_seqs = list(P1alignedtimefiltseqs_data.loc[:, 'Sequence_ids'])  
    P1Fitted_tfiltered_times = list(P1alignedtimefiltseqs_data.loc[:, 'Sequence_times'])    
    # fit by port1:
    Port1Pokes = port_fitted_poke_times(P1Fitted_tfiltered_seqs, P1Fitted_tfiltered_times, port1, 6)

    start_port = 0
    df = create_plotting_df(Port1Pokes, new_order)

    colors = ['#E6BFC4', '#285CA6', '#219CC1', '#A4D3B6', '#ECEFB6', '#E6BFC4']

    # plots:
    dx = "index"
    dy = "Time"
    ort = "h"
    pal = colors
    sigma = .2
    f, ax = plt.subplots(figsize=(20, 20))

    ax = sns.stripplot(x=dy, y=dx, data=df, palette=pal, edgecolor="white", size=5, jitter=0.4, zorder=0, orient=ort)

    ax.set_xlabel('Time from Initiation', fontsize=30)
    ax.set_ylabel('', fontsize=30)
    ax.tick_params(axis="y", labelsize=20)
    ax.tick_params(axis="x", labelsize=20)

    ax.set_xlim(0, 2)

    all_indexes = list(df.loc[:, 'index'])
    plottedports = []
    for x in all_indexes:
        if x not in plottedports:
            plottedports.append(x)

    lst = []
    colors = []
    for index, portname in enumerate(plottedports):
        if portname == start_port:
            lst.append('Origin')
            colors.append('k')
        elif portname == 5:
            lst.append('Other')
            colors.append('grey')
        else:
            lst.append('Port ' + str((portname + 1)))
            colors.append('red')

    ax.set_title('Origin = Port ' + str(start_port + 1), loc='left', fontsize=20, color='grey')

    ax.set_yticklabels(lst, fontsize=30)
    [t.set_color(i) for (i, t) in zip(colors, ax.yaxis.get_ticklabels())]

    SaveFig((CurrentAnimal + '_' + file + '_' + Experiment_type[0][2::] + '_' + 'Port1FittedStripPlot_PokeInTimes.png'), InputPathCurrent + '\\SingleSessions\\' + file + '\\')
    
    return df

def plot_transition_trajectories(timefiltseqs_data, Session_data, CurrentAnimal, file, Experiment_type, InputPathCurrent):
    # Pull data from dataframes
    TimeFiltered_seqs = list(timefiltseqs_data.loc[:, 'Sequence_ids'])
    port1 = list(Session_data.loc[:, 'Port1'])[0]
    port2 = list(Session_data.loc[:, 'Port2'])[0]
    port3 = list(Session_data.loc[:, 'Port3'])[0]
    port4 = list(Session_data.loc[:, 'Port4'])[0]
    port5 = list(Session_data.loc[:, 'Port5'])[0]

    # Set needed reference frame:
    TransitionTypesIndex = np.array([11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36, 37, 38, 41, 42, 43, 44, 45, 46, 47, 48, 51, 52, 53, 54, 55, 56, 57, 58, 61, 62, 63, 64, 65, 66, 67, 68, 71, 72, 73, 74, 75, 76, 77, 78, 81, 82, 83, 84, 85, 86, 87, 88])
    
    # Determine numbers of transitions (time filtered sequences + repeat pokes removed)
    trajects = []
    for inds, seqs in enumerate(TimeFiltered_seqs):
        seqs = literal_eval(seqs)  # convert back from weird df string conversion thing
        for ind, transits in enumerate(seqs):
            if not str(transits)[0] == str(transits)[1]:
                trajects = np.append(trajects, transits)
    transition_number = []
    for transit_types in TransitionTypesIndex:
        temp = (np.where(trajects == float(transit_types)))
        transition_number.append(len(temp[0]))
    
    # Set xy coords of each trajectory:
    x = ([1, 1], [1, 2], [1, 3], [1, 4], [1, 1], [1, 2], [1, 3], [1, 4],
         [2, 1], [2, 2], [2, 3], [2, 4], [2, 1], [2, 2], [2, 3], [2, 4],
         [3, 1], [3, 2], [3, 3], [3, 4], [3, 1], [3, 2], [3, 3], [3, 4],
         [4, 1], [4, 2], [4, 3], [4, 4], [4, 1], [4, 2], [4, 3], [4, 4],
         [1, 1], [1, 2], [1, 3], [1, 4], [1, 1], [1, 2], [1, 3], [1, 4],
         [2, 1], [2, 2], [2, 3], [2, 4], [2, 1], [2, 2], [2, 3], [2, 4],
         [3, 1], [3, 2], [3, 3], [3, 4], [3, 1], [3, 2], [3, 3], [3, 4],
         [4, 1], [4, 2], [4, 3], [4, 4], [4, 1], [4, 2], [4, 3], [4, 4])
    y = ([1, 0.5], [1, 1], [1, 1], [1, 1], [1, 2], [1, 2], [1, 2], [1, 2],
         [1, 1], [1, 0.5], [1, 1], [1, 1], [1, 2], [1, 2], [1, 2], [1, 2],
         [1, 1], [1, 1], [1, 0.5], [1, 1], [1, 2], [1, 2], [1, 2], [1, 2],
         [1, 1], [1, 1], [1, 1], [1, 0.5], [1, 2], [1, 2], [1, 2], [1, 2],
         [2, 1], [2, 1], [2, 1], [2, 1], [2, 2.5], [2, 2], [2, 2], [2, 2],
         [2, 1], [2, 1], [2, 1], [2, 1], [2, 2], [2, 2.5], [2, 2], [2, 2],
         [2, 1], [2, 1], [2, 1], [2, 1], [2, 2], [2, 2], [2, 2.5], [2, 2],
         [2, 1], [2, 1], [2, 1], [2, 1], [2, 2], [2, 2], [2, 2], [2, 2.5])

    # Generate plot:
    fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(1, 5, figsize=(75, 15))
    ax1.axis('off')
    ax2.axis('off')
    ax3.axis('off')
    ax4.axis('off')
    ax5.axis('off')
    ax1.text(0.7, 2.7, ('Transitions from port 1'), fontsize=25, color='k')
    ax2.text(0.7, 2.7, ('Transitions from port 2'), fontsize=25, color='k')
    ax3.text(0.7, 2.7, ('Transitions from port 3'), fontsize=25, color='k')
    ax4.text(0.7, 2.7, ('Transitions from port 4'), fontsize=25, color='k')
    ax5.text(0.7, 2.7, ('Transition from all ports'), fontsize=25, color='k')
    ax1.set_ylim([0, 5])
    ax1.set_xlim([0, 5])
    ax2.set_ylim([0, 5])
    ax2.set_xlim([0, 5])
    ax3.set_ylim([0, 5])
    ax3.set_xlim([0, 5])
    ax4.set_ylim([0, 5])
    ax4.set_xlim([0, 5])
    ax5.set_ylim([0, 5])
    ax5.set_xlim([0, 5])

    Colours = ['black', 'firebrick', 'grey']

    # Add port reference circles
    circles = [
        plt.Circle((1, 1), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((2, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((3, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((4, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((1, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((2, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((3, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((4, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((1, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((2, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((3, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((4, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((1, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((2, 2), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((3, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((4, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((1, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((2, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((3, 1), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((4, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((1, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((2, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((3, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((4, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((1, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((2, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((3, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((4, 1), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((1, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((2, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((3, 2), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((4, 2), 0.3, color=Colours[2], fill=False, linewidth=5),
        plt.Circle((1, 1), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((2, 1), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((3, 1), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((4, 1), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((1, 2), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((2, 2), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((3, 2), 0.3, color=Colours[0], fill=False, linewidth=5),
        plt.Circle((4, 2), 0.3, color=Colours[0], fill=False, linewidth=5)
    ]

    for i, circle in enumerate(circles):
        if i < 8:
            ax1.add_artist(circle)
        elif i < 16:
            ax2.add_artist(circle)
        elif i < 24:
            ax3.add_artist(circle)
        elif i < 32:
            ax4.add_artist(circle)
        else:
            ax5.add_artist(circle)

    # Plot lines
    for transition in range(0, 5):
        if transition == 0:
            port = port1
        if transition == 1:
            port = port2
        if transition == 2:
            port = port3
        if transition == 3:
            port = port4
        if port > 4:
            a = port - 4
            b = 2
        else:
            a = port
            b = 1
        REF_circle = plt.Circle((a, b), 0.3, color=Colours[1], fill=False, linewidth=5)

        if transition == 0:
            ax1.add_artist(REF_circle)
            normalise = max(transition_number[((port1 * 8) - 8):((port1 * 8) - 8) + 8])
            for i in range(((port1 * 8) - 8), ((port1 * 8) - 8) + 8):
                ax1.plot(x[i], y[i], markersize=130, markeredgewidth=0.001, color='firebrick', marker='None', linewidth=(5 * (transition_number[i] / normalise)))
                ax1.plot(x[i][1], y[i][1], 'o', color='firebrick', alpha=0.3, markersize=(60 * (transition_number[i] / normalise)))
        if transition == 1:
            ax2.add_artist(REF_circle)
            normalise = max(transition_number[((port2 * 8) - 8):((port2 * 8) - 8) + 8])
            for i in range(((port2 * 8) - 8), ((port2 * 8) - 8) + 8):
                ax2.plot(x[i], y[i], markersize=130, markeredgewidth=0.001, color='firebrick', marker='None', linewidth=(5 * (transition_number[i] / normalise)))
                ax2.plot(x[i][1], y[i][1], 'o', color='firebrick', alpha=0.3, markersize=(60 * (transition_number[i] / normalise)))
        if transition == 2:
            ax3.add_artist(REF_circle)
            normalise = max(transition_number[((port3 * 8) - 8):((port3 * 8) - 8) + 8])
            for i in range(((port3 * 8) - 8), ((port3 * 8) - 8) + 8):
                ax3.plot(x[i], y[i], markersize=130, markeredgewidth=0.001, color='firebrick', marker='None', linewidth=(5 * (transition_number[i] / normalise)))
                ax3.plot(x[i][1], y[i][1], 'o', color='firebrick', alpha=0.3, markersize=(60 * (transition_number[i] / normalise)))
        if transition == 3:
            ax4.add_artist(REF_circle)
            normalise = max(transition_number[((port4 * 8) - 8):((port4 * 8) - 8) + 8])
            for i in range(((port4 * 8) - 8), ((port4 * 8) - 8) + 8):
                ax4.plot(x[i], y[i], markersize=130, markeredgewidth=0.001, color='firebrick', marker='None', linewidth=(5 * (transition_number[i] / normalise)))
                ax4.plot(x[i][1], y[i][1], 'o', color='firebrick', alpha=0.3, markersize=(60 * (transition_number[i] / normalise)))

        if transition == 4:
            for i in range(len(transition_number)):
                normalise = max(transition_number)
                ax5.plot(x[i], y[i], markersize=130, markeredgewidth=0.001, color='firebrick', marker='None', linewidth=(5 * (transition_number[i] / normalise)))

    ax1.plot(1, 3.5, 'o', color='firebrick', alpha=0.3, markersize=60)
    ax1.text(1.3, 3.45, ('=  ' + str(normalise) + '  transitions'), fontsize=25, color='k')

    # Save
    SaveFig((CurrentAnimal + '_' + file + '_' + Experiment_type[0][2::] + '_' + 'TransitionVectors.png'), InputPathCurrent + '\\SingleSessions\\' + file + '\\')

def plot_transition_times(TransitionLatency_Tfilt, TransitionTypes_Tfilt, Transitions, CurrentAnimal, file, Experiment_type, InputPathCurrent):
    ## Plot 1
    #filter data by transitions and put into format that is matplotlib friendly
    Transition_Latencies = Find_Transition_times(TransitionLatency_Tfilt, TransitionTypes_Tfilt, Transitions)
    #plot:
    fig, [ax1, ax2, ax3, ax4] = plt.subplots(1, 4, figsize=(60, 15))
    bins = 60
    a = ax1.hist(Transition_Latencies[0], bins=bins, range=[0, 2], color='firebrick')
    b = ax2.hist(Transition_Latencies[1], bins=bins, range=[0, 2], alpha=1, color='grey')
    c = ax3.hist(Transition_Latencies[2], bins=bins, range=[0, 2], alpha=1, color='lightblue')
    d = ax4.hist(Transition_Latencies[3], bins=bins, range=[0, 2], alpha=1, color='pink')
    #set axis and labels 
    ax1.set_xlim([0, 2])
    ax1.set_ylim([0, (max(a[0]) + 5)])
    ax2.set_xlim([0, 2])
    ax2.set_ylim([0, (max(b[0]) + 5)])
    ax3.set_xlim([0, 2])
    ax3.set_ylim([0, (max(c[0]) + 5)])
    ax4.set_xlim([0, 2])
    ax4.set_ylim([0, (max(d[0]) + 5)])
    ax1.set_xlabel('Time', fontsize=25)
    ax1.set_ylabel('Transitions', fontsize=25)
    ax2.set_xlabel('Time', fontsize=25)
    ax2.set_ylabel('Transitions', fontsize=25)
    ax3.set_xlabel('Time', fontsize=25)
    ax3.set_ylabel('Transitions', fontsize=25)
    ax4.set_xlabel('Time', fontsize=25)
    ax4.set_ylabel('Transitions', fontsize=25)
    ax1.tick_params(axis="y", labelsize=20)
    ax2.tick_params(axis="y", labelsize=20)
    ax3.tick_params(axis="y", labelsize=20)
    ax4.tick_params(axis="y", labelsize=20)
    ax1.tick_params(axis="x", labelsize=20)
    ax2.tick_params(axis="x", labelsize=20)
    ax3.tick_params(axis="x", labelsize=20)
    ax4.tick_params(axis="x", labelsize=20)
    ax1.set_title(('Transition 1'), loc='left', fontsize=30, pad=30, color='k')
    ax2.set_title(('Transition 2'), loc='left', fontsize=30, pad=30, color='k')
    ax3.set_title(('Transition 3'), loc='left', fontsize=30, pad=30, color='k')
    ax4.set_title(('Transition 4'), loc='left', fontsize=30, pad=30, color='k')
    #median lines:
    a_median = np.median(Transition_Latencies[0])
    b_median = np.median(Transition_Latencies[1])
    c_median = np.median(Transition_Latencies[2])
    d_median = np.median(Transition_Latencies[3])
    ax1.plot([a_median, a_median], [(max(a[0]) + 5), 0], '--', linewidth=4, color='green')
    ax2.plot([b_median, b_median], [(max(b[0]) + 5), 0], '--', linewidth=4, color='green')
    ax3.plot([c_median, c_median], [(max(c[0]) + 5), 0], '--', linewidth=4, color='green')
    ax4.plot([d_median, d_median], [(max(d[0]) + 5), 0], '--', linewidth=4, color='green')
    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
    #save
    SaveFig((CurrentAnimal + '_' + file + '_' + Experiment_type[0][2::] + '_' + 'TransitionTimes1_PokeOut_PokeIn.png'), InputPathCurrent + '\\SingleSessions\\' + file + '\\')

            

def plot_transition_times2(TransitionLatency_Tfilt, TransitionTypes_Tfilt, Transitions, CurrentAnimal, file, Experiment_type, InputPathCurrent):
    # pull out data in df friendly format and save into dataframe
    Transition_Latencies, transition_type = Find_Transition_times_sns(TransitionLatency_Tfilt, TransitionTypes_Tfilt, Transitions)
    if transition_type:  # if it's not empty
        df = pd.DataFrame({'Transition': transition_type, 'Time': Transition_Latencies})
        # plots:
        dx = "Transition"
        dy = "Time"
        ort = "v"
        pal = "Set2"
        sigma = .2
        f, ax = plt.subplots(figsize=(30, 15))
        ax = pt.half_violinplot(x=dx, y=dy, data=df, palette=pal, bw=.2, cut=0., scale="area", width=.6, inner=None, orient=ort)
        ax = sns.stripplot(x=dx, y=dy, data=df, palette=pal, edgecolor="white", size=3, jitter=1, zorder=0, orient=ort)
        ax = sns.boxplot(x=dx, y=dy, data=df, color="grey", width=.15, zorder=10, showcaps=True, boxprops={'facecolor': 'none', "zorder": 10}, showfliers=False, whiskerprops={'linewidth': 2, "zorder": 10}, saturation=1, orient=ort)
        # set labels
        ax.set_ylabel('Time', fontsize=30)
        ax.set_xlabel('', fontsize=30)
        ax.tick_params(axis="y", labelsize=15)

        splitdata = dict(tuple(df.groupby('Transition')))
        labels = []
        for i in splitdata:
            labels.append('Transition ' + str(i))
        ax.set_xticklabels(labels, fontsize=20)
        # save
        SaveFig((CurrentAnimal + '_' + file + '_' + Experiment_type[0][2::] + '_' + 'TransitionTimes2_PokeOut_PokeIn.png'), InputPathCurrent + '\\SingleSessions\\' + file + '\\')
