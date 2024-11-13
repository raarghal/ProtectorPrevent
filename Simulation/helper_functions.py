import numpy as np

def calc_r_0(pop_sizes, C, state_lengths, transmission_rates, trans_probabilities):
    y = []
    for cluster in range(C.shape[0]):
        symp_length = state_lengths["P"][cluster]+state_lengths["I"][cluster]
        avg_inf_length = trans_probabilities["symp"][cluster]*symp_length + (1 - trans_probabilities["symp"][cluster])*state_lengths["A"][cluster]
        y.append(avg_inf_length)
    next_gen_mat = np.diag(transmission_rates) @ (C @ np.diag(y))
    return max(np.absolute(np.linalg.eigvals(next_gen_mat)))

def scale_contacts(pop_sizes, N, contact_matrix, state_lengths, transmission_rates, trans_probabilities, npi_efficacy, variant, r_0 = None):
    C = contact_matrix
#    C = (C + C.T)/2
    #for i in range(C.shape[0]):
#        for j in range(C.shape[0]):
#            C[i,j] = C[i,j]/pop_sizes[j]

    if (r_0 is not None):
        #print("hi")
        #print(calc_r_0(pop_sizes, C, state_lengths, transmission_rates, trans_probabilities))
        #print(C)
        C *= r_0/(calc_r_0(pop_sizes, C, state_lengths, transmission_rates, trans_probabilities))
        #print(C)
        #print(calc_r_0(pop_sizes, C, state_lengths, transmission_rates, trans_probabilities))
    C = np.diag([1 - el for el in npi_efficacy]) @ C
    #C = np.diag([1/N]*3) @ C
    return C

def epi_step(pop_sizes, current_state, C, state_lengths, trans_probabilities, transmission_rates, variant, vax_efficacy, timestep, r_0):
    new_infs = [0]*current_state.shape[0]
    new_hosps = [0]*current_state.shape[0]
    next_state = current_state
    for cluster_1 in range(current_state.shape[0]):
        num_exposures = 0.0
        num_vax_exposures = 0.0
        for cluster_2 in range(current_state.shape[0]):
            for infectious_state in range(2,5,1):
                contacts = C[cluster_1, cluster_2] * current_state[cluster_1, 0] * current_state[cluster_2, infectious_state]/pop_sizes[cluster_2]
                #print(contacts)
                exposures = transmission_rates[cluster_1] * contacts

                #print(num_exposures)
                vax_contacts = C[cluster_1, cluster_2] * current_state[cluster_1, 7] * current_state[cluster_2, infectious_state]/pop_sizes[cluster_2]
                vax_exposures = transmission_rates[cluster_1] * vax_contacts

                num_exposures += exposures * variant["exp"][cluster_1]
                num_vax_exposures += vax_exposures * variant["exp"][cluster_1] * (1-vax_efficacy[cluster_1])
        #print(num_exposures * timestep)
        new_infs[cluster_1] = timestep*(num_exposures + num_vax_exposures)


        next_state[cluster_1, 0] -= num_exposures * timestep
        next_state[cluster_1, 1] += num_exposures * timestep


        num_E_to_P = trans_probabilities["symp"][cluster_1] * variant["symp"][cluster_1] * (1/state_lengths["E"][cluster_1]) * current_state[cluster_1, 1]
        next_state[cluster_1, 1] -= num_E_to_P * timestep
        next_state[cluster_1, 2] += num_E_to_P * timestep

        num_E_to_A = (1 - trans_probabilities["symp"][cluster_1] * variant["symp"][cluster_1]) * (1/state_lengths["E"][cluster_1]) * current_state[cluster_1, 1]
        next_state[cluster_1, 1] -= num_E_to_A * timestep
        next_state[cluster_1, 3] += num_E_to_A * timestep

        num_P_to_I = current_state[cluster_1, 2] * (1/state_lengths["P"][cluster_1])
        next_state[cluster_1, 2] -= num_P_to_I * timestep
        next_state[cluster_1, 4] += num_P_to_I * timestep

        num_A_to_R = current_state[cluster_1, 3] * (1/state_lengths["A"][cluster_1])
        next_state[cluster_1, 3] -= num_A_to_R * timestep
        next_state[cluster_1, 8] += num_A_to_R * timestep

        num_I_to_L = current_state[cluster_1, 4] * (1/state_lengths["I"][cluster_1])
        next_state[cluster_1, 4] -= num_I_to_L * timestep
        next_state[cluster_1, 5] += num_I_to_L * timestep

        num_L_to_H = current_state[cluster_1, 5] * (1/state_lengths["L"][cluster_1]) * trans_probabilities["hosp"][cluster_1] * variant["hosp"][cluster_1]
        new_hosps[cluster_1] = timestep*(num_L_to_H)
        next_state[cluster_1, 5] -= num_L_to_H * timestep
        next_state[cluster_1, 6] += num_L_to_H * timestep

        num_L_to_R = current_state[cluster_1, 5] * (1/state_lengths["L"][cluster_1]) * (1 - trans_probabilities["hosp"][cluster_1] * variant["hosp"][cluster_1])
        next_state[cluster_1, 5] -= num_L_to_R * timestep
        next_state[cluster_1, 8] += num_L_to_R * timestep

        num_H_to_R = current_state[cluster_1, 6] * (1/state_lengths["H"][cluster_1]) * (1 - trans_probabilities["dec"][cluster_1] * variant["dec"][cluster_1])
        next_state[cluster_1, 6] -= num_H_to_R * timestep
        next_state[cluster_1, 8] += num_H_to_R * timestep

        num_H_to_D = current_state[cluster_1, 6] * (1/state_lengths["H"][cluster_1]) * trans_probabilities["dec"][cluster_1] * variant["dec"][cluster_1]
        next_state[cluster_1, 6] -= num_H_to_D * timestep
        next_state[cluster_1, 9] += num_H_to_D * timestep

        next_state[cluster_1, 7] -= num_vax_exposures * timestep
        next_state[cluster_1, 1] += num_vax_exposures * timestep
        #print(next_state)
        #print(sum(sum(next_state)))
    return next_state, new_hosps, new_infs

def vax_step(N, current_state, vax_rate, vax_efficacy, strategy, timestep):
    next_state = current_state
    vax_remain = vax_rate * N * timestep
    switcher = {
        "high risk": [1,2,0],
        "high contact": [2,1,0],
        "uniform": [0,1,2],
        "hchrc": [3,2,1,0],
        "hchrr": [3,1,2,0],
        "chchrr": [2,3,1,0],
        "rhchrc": [1,3,2,0],
        "rchchr": [1,2,3,0],
        "crhchr": [2,1,3,0],
        "none": [0,1,2]
    }
    #print(vax_remain)
    order = switcher[strategy]
    rates = np.array([[0.0]*4])
    nums = np.array([[0.0]*4])
    if(strategy != "uniform" and strategy !="none"):
        for cluster in order:
            if (current_state[cluster,0] >= vax_remain):
                next_state[cluster,0] -= vax_remain
                next_state[cluster,7] += vax_remain
                rates[0,cluster] = vax_remain / current_state[cluster,0]
                nums[0,cluster] = vax_remain
                vax_remain = 0

                break
            else:
                next_state[cluster,0] = 0
                next_state[cluster,7] += current_state[cluster,0]
                rates[0,cluster] = 1.0
                nums[0,cluster] = current_state[cluster,0]
                vax_remain -= current_state[cluster,0]
    elif(strategy == "none"):
        rates = [0]*4
        nums = [0]*4
    else:
        if (current_state[0,0]+current_state[1,0]+current_state[2,0]>0):
            rate = vax_remain/(current_state[0,0]+current_state[1,0]+current_state[2,0])
            rates[0,:] = [rate]*4
            np.multiply(rates[0,:],current_state[:,0],nums[0,:])
            next_state[:,0] -= nums[0,:]
            next_state[:,7] += nums[0,:]
        else:
            rates = [0]*4
            nums = [0]*4
    #print(next_state)
    return next_state, rates, nums
