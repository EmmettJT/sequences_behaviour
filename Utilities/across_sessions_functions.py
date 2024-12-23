import scipy.io
# import os, importlib
import matplotlib.pyplot as plt
import statistics
import scipy.stats
import matplotlib.patches as mpatches
import seaborn as sns
import os
import tkinter as tk
from matplotlib import gridspec
import matplotlib.colors
import numpy as np
import pandas as pd
# import ptitprince as pt
from ast import literal_eval
from matplotlib import cm
# from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pickle
import shutil
from IPython.display import HTML

plt.style.use('default')

def convolve_movmean(y,N):
    y_padded = np.pad(y, (N//2, N-1-N//2), mode='edge')
    y_smooth = np.convolve(y_padded, np.ones((N,))/N, mode='valid') 
    return y_smooth

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(matplotlib.colors.to_rgb(c1))
    c2=np.array(matplotlib.colors.to_rgb(c2))
    return matplotlib.colors.to_hex((1-mix)*c1 + mix*c2)

def SaveFig(file_name,figure_dir):
    if not os.path.isdir(figure_dir):
        os.makedirs(figure_dir)
    plt.savefig(figure_dir + file_name, bbox_inches='tight')
    plt.close()
    
def SaveFig_hold(file_name,figure_dir):
    if not os.path.isdir(figure_dir):
        os.makedirs(figure_dir)
    plt.savefig(figure_dir + file_name, bbox_inches='tight')

def DetermineTransitionScores(TransitionTypes_Tfilt,current_LEDi,current_Rval,Transitions,port):
    transition_scores = []
    transition_LEDi = []
    transition_Rval = []
    for index,transit in enumerate(TransitionTypes_Tfilt):
        if str(transit)[0] == str(Transitions[port-1])[0]:
            if transit == Transitions[port-1]:
                transition_scores = transition_scores + [1]
                transition_LEDi = transition_LEDi + [literal_eval(current_LEDi[index])]
                transition_Rval = transition_Rval + [literal_eval(current_Rval[index])]

            else:
                transition_scores = transition_scores+ [0]
                transition_LEDi = transition_LEDi + [literal_eval(current_LEDi[index])]
                transition_Rval = transition_Rval + [literal_eval(current_Rval[index])]
    return transition_scores, transition_LEDi, transition_Rval

def convolve_movmean(y,N):
    y_padded = np.pad(y, (N//2, N-1-N//2), mode='edge')
    y_smooth = np.convolve(y_padded, np.ones((N,))/N, mode='valid') 
    return y_smooth

def Find_Transition_times(TransitionLatency_Tfilt,TransitionTypes_Tfilt,Transitions):
    transit_times = []
    for Transits in Transitions:
        transit_times_temp = []
        for ind, transition_pair in enumerate(TransitionTypes_Tfilt):
            if transition_pair == Transits:
                transit_times_temp = transit_times_temp + [TransitionLatency_Tfilt[ind]]
                
        transit_times = transit_times + [transit_times_temp]   
    return(transit_times)

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

def mklist(n):
    for _ in range(n):
        yield []
        
#### align to first port pokes and remove single transitions (these dont count as sequences)

def aligntofirstpokeandremovesingletransits(timesplitseqs,timesplitlatencies):
    
    newseqs = []
    newlatencies = []
    # align to first poke:
    for index_1,fragments in enumerate(timesplitseqs):
        current_newseqs = []
        current_newlatencies = []
        count = -1
        seqs = False
        for index_2,sequence in enumerate(fragments):
            for index_3,transit in enumerate(sequence):
                if not str(transit)[0] == str(transit)[1]: # remove repeat pokes
                    if str(transit)[0] == '2':
                        seqs = True
                        current_newseqs = current_newseqs + [[]]
                        current_newlatencies = current_newlatencies + [[]]
                        count = count + 1
                        current_newseqs[count] = current_newseqs[count] + [transit]
                        current_newlatencies[count] = current_newlatencies[count] + [timesplitlatencies[index_1][index_2][index_3]]
                    elif seqs == True:
                        current_newseqs[count] = current_newseqs[count] + [transit]   
                        current_newlatencies[count] = current_newlatencies[count] + [timesplitlatencies[index_1][index_2][index_3]]
            seqs = False
 
        newseqs = newseqs + [current_newseqs]
        newlatencies = newlatencies + [current_newlatencies]
    return(newseqs,newlatencies)

def sequence_contains_sequence(haystack_seq, needle_seq):
    for i in range(0, len(haystack_seq) - len(needle_seq) + 1):
        if needle_seq == haystack_seq[i:i+len(needle_seq)]:
            return True
    return False
            
def parts(list_, indices):
    indices = [0]+indices+[len(list_)]
    return [list_[v:indices[k+1]] for k, v in enumerate(indices[:-1])]

def RemoveSlowSequences(split,split2):
    timefiltered_split = []
    for i,item in enumerate(split2):
        if item[0] == 1:
            timefiltered_split = timefiltered_split + [split[i]]

    return(timefiltered_split)

def determine_transition_score(fragments,start_port,end_port):
    correct = []
    for fragment in fragments:
        for transition in fragment:
            if str(transition)[0] == str(start_port):
                if str(transition)[1] == str(end_port):
                    correct = correct + [1]
                else:
                    correct = correct + [0]
    return(correct)

def generate_processed_transitiontimesdataframe(processed_seqs,processed_latencies,counter):

    count = counter
    transits= []
    trial_number= []
    for fragment in processed_seqs:
        count = count + 1
        if len(fragment) > 0:
            for sequence in fragment:
                for transit in sequence:
                    trial_number = trial_number + [count]
                    transits = transits + [transit]
        else: ### deals with cases where there are no good transitions in a trial 
            transits = transits + ['nan']
            trial_number = trial_number + [count]

    times = []
    for fragment in processed_latencies:
        if len(fragment) > 0:
            for sequence in fragment:
                for time in sequence:
                    times = times + [time]
        else:
            times = times + ['nan']

    Processesed_Transition_Latencies = pd.DataFrame({'Trial': trial_number, 'Transitions' : transits,'Latencies' : times})

    return(Processesed_Transition_Latencies,count)

def findTransitionTimes(Transition_id,all_sessions_latency_data):
    Transition = np.array(pd.concat(all_sessions_latency_data)['Latencies'])[np.where(np.array(pd.concat(all_sessions_latency_data)['Transitions'])==Transition_id)]
    Trials = np.array(pd.concat(all_sessions_latency_data)['Trial'])[np.where(np.array(pd.concat(all_sessions_latency_data)['Transitions'])==Transition_id)]
    return(Transition,Trials)

def find_mean_transiton_times_by_trail(transitions,trials):
    Mean_Transitions_bytrial = []
    ref_trials = []
    for i in range(len(transitions)):
        df =  pd.DataFrame({'Trials': trials[i], 'Transitions' : transitions[i]})
        split_by_trials = dict(tuple(df.groupby('Trials')))

        c_Mean_Transitions_bytrial = []
        c_ref_trials= []
        for item in split_by_trials:
            current_mean = np.mean(np.array(split_by_trials[item]['Transitions']))
            c_Mean_Transitions_bytrial = c_Mean_Transitions_bytrial + [current_mean]
            c_ref_trials = c_ref_trials + [item]

        Mean_Transitions_bytrial = Mean_Transitions_bytrial + [c_Mean_Transitions_bytrial]
        ref_trials = ref_trials + [c_ref_trials]

    return(Mean_Transitions_bytrial,ref_trials)