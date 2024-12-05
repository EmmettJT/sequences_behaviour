# sequences_behaviour
scripts for running and post-hoc analysis of the 8 poke sequence behavioural task 


# Analysis of behavioural data for the Seqeuence task 

See our publication for more details: [Replay of Procedural Experience is Independent of the Hippocampus](https://www.biorxiv.org/content/10.1101/2024.06.05.597547v1.full.pdf).

![Task Schematic](images/schematic.png)

## Overview

This repository includes:
- A version of the task in BPOD file format (Sequence_Automated)
- Tidied* scripts for processing BPOD output
- Some example scripts for analysis of tracking and histology data 

* These scripts should be useful and readable but are not fully polished.

## Getting Started

Ensure you have the following software installed:
- [Git](https://git-scm.com/)
- [Python](https://www.python.org/downloads/)  (Version used: 3.12.3)
- Jupyter notebook
- Necessary Python libraries: this pipeline worked with the environment I have listed in requirements.txt (though some of these packages may be redundant)  

### Data Pipeline

**1_process_bpod_output**
Notebook script which takes the matlab output from running the behaviour and produces more readable processed data 

**2_analyse_individual_sessions**
Takes output files from process_bpod_output and returns analysis plots for each session

**3_analyse_across_sessions**
same as step 2 but produces summary data across the sessions for each animal 

**AcrossAnimal_analysis**
takes the outputs from the first 3 steps and creates plots which compare the task abilities of multiple animals 
