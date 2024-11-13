from datetime import datetime
import numpy as np
import pandas as pd
import xlsxwriter as xl
import matplotlib.pyplot as plt
from operator import add
from openpyxl import load_workbook
from itertools import islice

from helper_functions import *

def runSim(init_state, contact_matrix, state_lengths, trans_probabilities, transmission_rates, vax_rate, vax_efficacy, npi_efficacy, variant, strategy, r_0 = None):
    N = sum(sum(init_state))
    pop_sizes = np.sum(init_state,axis=1)
    #print(pop_sizes)
    C = scale_contacts(pop_sizes, N, contact_matrix, state_lengths, transmission_rates, trans_probabilities, npi_efficacy, variant, r_0)
    print(C)

    T = 600.0
    timestep = 0.01

    num_clusters = init_state.shape[0]
    num_states = init_state.shape[1]

    sim_results = np.tile(0.0,(round(T/timestep)+1, num_clusters, num_states))
    sim_results[0,:,:] = init_state
    vax_rates = np.tile(0.0,(round(T/timestep)+1, num_clusters))
    vax_nums = np.tile(0.0,(round(T/timestep)+1, num_clusters))
    total_infs = [0]*num_clusters
    total_hosps = [0]*num_clusters
    cumul_infs = 0
    cumul_hosps = 0

    for i in range(1,round(T/timestep)+1,1):
        #print(C)

        sim_results[i,:,:], new_hosps, new_infs = epi_step(pop_sizes, sim_results[i-1], C, state_lengths, trans_probabilities, transmission_rates, variant, vax_efficacy, timestep, r_0)
        total_infs = [sum(x) for x in zip(total_infs, new_infs)]
        total_hosps = [sum(x) for x in zip(total_hosps, new_hosps)]
        cumul_infs += sum(sim_results[i,:,4])*timestep
        cumul_hosps += sum(sim_results[i,:,6])*timestep
        #print(sim_results[i,:,:])
        #print(abs(sum(sum(sim_results[i,:,:]))))
        #sim_results[sim_results < 0] = 0
        #print(abs(sum(sum(sim_results[i,:,:]))))
        sim_results[i,:,:], vax_rates[i,:], vax_nums[i,:] = vax_step(N, sim_results[i], vax_rate, vax_efficacy, strategy, timestep)
        #print(abs(sum(sum(sim_results[i,:,:]))))
        #print(sim_results[i,:,:])
        #sim_results[sim_results < 0] = 0
        #print(abs(sum(sum(sim_results[i,:,:]))))
        #print(sim_results[i,:,:])
        #if(strategy == "uniform"):
            #print(sim_results[i,:,9])
        #if (abs(sum(sum(sim_results[i,:,:]))-N)>10):
            #print(sum(sum(sim_results[i,:,:])))
            #raise Exception("something is wrong")
    return sim_results, vax_rates, vax_nums, total_hosps, total_infs,cumul_infs,cumul_hosps
