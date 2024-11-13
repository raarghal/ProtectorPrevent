from datetime import datetime
import numpy as np
import pandas as pd
import xlsxwriter as xl
import matplotlib.pyplot as plt
from operator import add
from openpyxl import load_workbook
from itertools import islice
from itertools import cycle
import seaborn as sns
from RunSim import *
from helper_functions import *

vax_efficacy_trials = [0.9, 0.8, 0.7]
npi_efficacy_trials = [0, 0.3, 0.6]
initial_infections_trials = [15]#[2,5,15,20]#[0.05, 0.1, 0.5, 1]#[0.1, 0.5, 1, 5]#[1,5,10,20]
variant_trials = [1]#[1,1.6,2.2, 2.5]
r_0_vals = [2]#[1, 1.5,2,3]
HC_pop = [0.1,0.15,0.2]
settings = ['USA']#['USA']#['USA', 'LMIC', 'Prison', 'Retirement Home']


for setting in settings:
    results = xl.Workbook('Data/fullSim_results_'+setting+'_{}.xlsx'.format(datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")))
    death_reduction = results.add_worksheet('Death Reduction')
    death_reduction.write(0,0,'Vaccine Efficacy')
    death_reduction.write(0,1,'NPI Efficacy')
    death_reduction.write(0,2,'Initial Infections')
    death_reduction.write(0,3,'Variant Parameter')
    death_reduction.write(0,4, 'R_0')
    death_reduction.write(0,5, 'Setting')
    death_reduction.write(0,6, 'High Contact Population')
    death_reduction.write(0,7,'Uniform Vaccination Deaths')
    death_reduction.write(0,8,'Uniform Vaccination Hospitalizations')
    death_reduction.write(0,9,'Uniform Vaccination Infections')
    death_reduction.write(0,10,'Uniform Cumulative Hospitalizations')
    death_reduction.write(0,11,'Uniform Cumulative Infections')
    death_reduction.write(0,12,'High Risk Vaccination Deaths')
    death_reduction.write(0,13,'High Risk Vaccination Hospitalizations')
    death_reduction.write(0,14,'High Risk Vaccination Infections')
    death_reduction.write(0,15,'High Risk Cumulative Hospitalizations')
    death_reduction.write(0,16,'High Risk Cumulative Infections')
    death_reduction.write(0,17,'High Contact Vaccination Deaths')
    death_reduction.write(0,18,'High Contact Vaccination Hospitalizations')
    death_reduction.write(0,19,'High Contact Vaccination Infections')
    death_reduction.write(0,20,'High Contact Cumulative Hospitalizations')
    death_reduction.write(0,21,'High Contact Cumulative Infections')
    death_reduction.write(0,22,'No Vaccination Baseline Deaths')
    death_reduction.write(0,23,'No Vaccination Baseline Hospitalizations')
    death_reduction.write(0,24,'No Vaccination Baseline Infections')
    death_reduction.write(0,25,'No Vaccination Cumulative Hospitalizations')
    death_reduction.write(0,26,'No Vaccination Cumulative Infections')
    record = 1
    for vax_effi in vax_efficacy_trials:
        for npi_effi in npi_efficacy_trials:
            for init_inf in initial_infections_trials:
                for var_param in variant_trials:
                    for r in r_0_vals:
                        for f_hc in HC_pop:
                            print("Running Trial " + str(record) + " ..." )
                            death_reduction.write(record,0,vax_effi)
                            death_reduction.write(record,1,npi_effi)
                            death_reduction.write(record,2,init_inf)
                            death_reduction.write(record,3,var_param)
                            death_reduction.write(record,4,r)
                            death_reduction.write(record,5,setting)
                            death_reduction.write(record,6,f_hc)

                            if (setting == 'USA'):
                                f_high_risk = 0.168
                                f_high_contact = f_hc
                                contact_matrix = np.array([
                                    [0.165, 0.1, 0.175],
                                    [0.1, 0, 0.002],
                                    [0.175, 0.002, 0.132]
                                ])
                                vax_rate = 0.005

                            elif (setting == 'FL'):
                                f_high_risk = 0.201
                                f_high_contact = 0.4/2
                                contact_matrix = np.array([
                                    [0.165, 0.1, 0.175],
                                    [0.1, 0, 0.002],
                                    [0.175, 0.002, 0.132]
                                ])
                                vax_rate = 0.005

                            elif (setting == 'AZ'):
                                f_high_risk = 0.171
                                f_high_contact = 0.393/2
                                contact_matrix = np.array([
                                    [0.165, 0.1, 0.175],
                                    [0.1, 0, 0.002],
                                    [0.175, 0.002, 0.132]
                                ])
                                vax_rate = 0.005

                            elif (setting == 'UT'):
                                f_high_risk = 0.108
                                f_high_contact = 0.446/2
                                contact_matrix = np.array([
                                    [0.165, 0.1, 0.175],
                                    [0.1, 0, 0.002],
                                    [0.175, 0.002, 0.132]
                                ])
                                vax_rate = 0.005

                            elif (setting == 'ND'):
                                f_high_risk = 0.149
                                f_high_contact = 0.572/2
                                contact_matrix = np.array([
                                    [0.165, 0.1, 0.175],
                                    [0.1, 0, 0.002],
                                    [0.175, 0.002, 0.132]
                                ])
                                vax_rate = 0.005

                            elif (setting == 'LMIC'):
                                f_high_risk = 0.10
                                f_high_contact = 0.20
                                contact_matrix = np.array([
                                    [0.165, 0.029, 0.175],
                                    [0.029, 0, 0.01],
                                    [0.175, 0.01, 0.132]
                                ])
                                vax_rate = 0.002

                            elif (setting == 'Prison'):
                                f_high_risk = 0.08
                                f_high_contact = 0.18
                                contact_matrix = np.array([
                                    [0.165, 0.165, 0.175],
                                    [0.165, 0.165, 0.175],
                                    [0.175, 0.175, 0.175]
                                ])
                                vax_rate = 0.005
                                var_param = 1.8

                            elif (setting == 'Retirement Home'):
                                f_high_risk = 0.60
                                f_high_contact = 0.20
                                contact_matrix = np.array([
                                    [0.165, 0.165, 0.175],
                                    [0.165, 0.002, 0.175],
                                    [0.175, 0.175, 0.132]
                                ])
                                vax_rate = 0.005


                            f_base = 1-f_high_risk-f_high_contact
                            N = 2000

                            init_state = np.array([
                                [f_base*(N-init_inf), f_base*init_inf, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                [f_high_risk*(N-init_inf), f_high_risk*init_inf, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                [f_high_contact*(N-init_inf), f_high_contact*init_inf, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                            ])

                            pop_sizes = np.sum(init_state,axis=1)

                            #contact_matrix = np.array([
                            #    [5.0, 1.5, 15.0/6.0],
                            #    [3.0, 2.0, 2.0/3.0],
                            #    [15.0, 2.0, 5.0]
                            #])

                            #contact_matrix = np.array([
                            #    [0.12, 0.09, 0.264],
                            #    [0.09, 0.09, 0.09],
                            #    [0.264, 0.110, 0.572]
                            #])



                            #contact_matrix *= 2

                            state_lengths = {
                                "E": [4.0, 4.0, 4.0],
                                "P": [2.0, 2.0, 2.0],
                                "A": [10.0, 10.0, 10.0],
                                "I": [3.0, 3.0, 3.0],
                                "L": [3.0, 3.0, 3.0],
                                "H": [11.0, 11.0, 11.0],
                            }

                            trans_probabilities = {
                                "symp": [0.4, 0.8, 0.4],
                                "hosp": [0.1, 0.3, 0.1],
                                "dec": [0.01, 0.1, 0.01],
                            }

                            #trans_probabilities = {
                            #    "symp": [0.4, 0.8, 0.4],
                            #    "hosp": [0.25, 0.375, 0.25],
                            #    "dec": [0.1, 0.33, 0.1],
                            #}

                            transmission_rates = np.array([
                                0.4, 0.8, 0.4
                            ])



                            vax_efficacy = np.array([
                                vax_effi, vax_effi, vax_effi
                            ])

                            npi_efficacy = np.array([
                                npi_effi, npi_effi, npi_effi
                            ])

                            variant = {
                                "exp": [var_param, var_param, var_param],
                                "symp": [var_param, var_param, var_param],
                                "hosp": [var_param, var_param, var_param],
                                "dec": [var_param, var_param, var_param]
                            }

                            baseline_results, baseline_rates, baseline_nums, baseline_hosps, baseline_infs, bcumul_infs,bcumul_hosps = runSim(init_state, contact_matrix, state_lengths, trans_probabilities, transmission_rates, vax_rate, vax_efficacy, npi_efficacy, variant, "none", r)
                            uniform_results, uniform_rates, uniform_nums, uniform_hosps, uniform_infs, ucumul_infs,ucumul_hosps = runSim(init_state, contact_matrix, state_lengths, trans_probabilities, transmission_rates, vax_rate, vax_efficacy, npi_efficacy, variant, "uniform", r)
                            high_risk_results, high_risk_rates, high_risk_nums, high_risk_hosps, high_risk_infs, rcumul_infs,rcumul_hosps = runSim(init_state, contact_matrix, state_lengths, trans_probabilities, transmission_rates, vax_rate, vax_efficacy, npi_efficacy, variant, "high risk", r)
                            high_contact_results, high_contact_rates, high_contact_nums, high_contact_hosps, high_contact_infs, ccumul_infs,ccumul_hosps = runSim(init_state, contact_matrix, state_lengths, trans_probabilities, transmission_rates, vax_rate, vax_efficacy, npi_efficacy, variant, "high contact", r)
                            deaths = np.array([uniform_results[-1,:,9].sum(), high_risk_results[-1,:,9].sum(), high_contact_results[-1,:,9].sum(), baseline_results[-1,:,9].sum()])
                            death_reduction.write(record,7,deaths[0])
                            death_reduction.write(record,8,sum(uniform_hosps))
                            death_reduction.write(record,9,sum(uniform_infs))
                            death_reduction.write(record,10,ucumul_hosps)
                            death_reduction.write(record,11,ucumul_infs)
                            death_reduction.write(record,12,deaths[1])
                            death_reduction.write(record,13,sum(high_risk_hosps))
                            death_reduction.write(record,14,sum(high_risk_infs))
                            death_reduction.write(record,15,rcumul_hosps)
                            death_reduction.write(record,16,rcumul_infs)
                            death_reduction.write(record,17,deaths[2])
                            death_reduction.write(record,18,sum(high_contact_hosps))
                            death_reduction.write(record,19,sum(high_contact_infs))
                            death_reduction.write(record,20,ccumul_hosps)
                            death_reduction.write(record,21,ccumul_infs)
                            death_reduction.write(record,22,deaths[3])
                            death_reduction.write(record,23,sum(baseline_hosps))
                            death_reduction.write(record,24,sum(baseline_infs))
                            death_reduction.write(record,25,bcumul_hosps)
                            death_reduction.write(record,26,bcumul_infs)
                            #death_rates = deaths
                            #print(deaths)
                            #np.divide(deaths, sum(pop_sizes), death_rates)
                            percent_reduc = -100*(deaths - deaths[3])/deaths[3]
                            record += 1
    #print(uniform_rates)
    results.close()
"""
objects = ('Uniform', 'High Risk', 'High Contact', 'None')
y_pos = np.arange(len(objects))
plt.bar(y_pos, percent_reduc, align = 'center', alpha = 0.5)
plt.xticks(y_pos, objects)
plt.ylabel('% Reduction in Deaths')
plt.title("Mortality")
plt.show()

t = range(0,60001,1)
plt.plot(t, baseline_results[:,0,0], label='S')
plt.plot(t, baseline_results[:,0,1], label='E')
plt.plot(t, baseline_results[:,0,2], label='P')
plt.plot(t, baseline_results[:,0,3], label='A')
plt.plot(t, baseline_results[:,0,4], label='I')
plt.plot(t, baseline_results[:,0,5], label='L')
plt.plot(t, baseline_results[:,0,6], label='H')
plt.plot(t, baseline_results[:,0,7], label='V')
plt.plot(t, baseline_results[:,0,8], label='R')
plt.plot(t, baseline_results[:,0,9], label='D')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('# of People')
plt.title('Epidemic Process w/ No Vaccination')
plt.show()
#print(uniform_results[3000,:,:])
#print(np.sum(uniform_results[3000,:,:],axis=1))

#t = range(0,10001,1)
#plt.plot(t, high_risk_nums[:,0], label='Low contact, low risk')
#plt.plot(t, high_risk_nums[:,1], label='low contact, high risk')
#plt.plot(t, high_risk_nums[:,2], label='high contact')
#plt.plot(t, np.add(uniform_nums[:,0],np.add(uniform_nums[:,1],uniform_nums[:,2])), label='total')

#plt.legend(loc='best')
#plt.xlabel('Time')
#plt.ylabel('High Risk Vax')
#plt.title('Epidemic')
#plt.show()
"""
uniform_deaths = uniform_results[-1,:,9]
high_risk_deaths = high_risk_results[-1,:,9]
high_contact_deaths = high_contact_results[-1,:,9]
baseline_deaths = baseline_results[-1,:,9]
#print(deaths)
hosps = [sum(uniform_hosps), sum(high_risk_hosps), sum(high_contact_hosps), sum(baseline_hosps)]
infs = [sum(uniform_infs), sum(high_risk_infs), sum(high_contact_infs), sum(baseline_infs)]
years_lost = [37,7,47]
YLLs = [sum([x*y for x,y in zip(years_lost,uniform_deaths)]), sum([x*y for x,y in zip(years_lost,high_risk_deaths)]), sum([x*y for x,y in zip(years_lost,high_contact_deaths)]), sum([x*y for x,y in zip(years_lost,baseline_deaths)])]
#print(uniform_results[3000,:,:])
percent_reduc_deaths = -100*(deaths - deaths[3])/deaths[3]
#print(percent_reduc_deaths[0:3])
percent_reduc_hosps = -100*(hosps - hosps[3])/hosps[3]
percent_reduc_infs = -100*(infs - infs[3])/infs[3]
percent_reduc_YLL = -100*(YLLs - YLLs[3])/YLLs[3]
#print("uniform deaths:")
#print(uniform_results[-1,:,9])
#print("high risk deaths:")
#print(high_risk_results[-1,:,9])
#print("high contact deaths")
#print(high_contact_results[-1,:,9])
#print("baseline deaths")
#print(baseline_results[-1,:,9])
redux = pd.DataFrame()
strats = cycle(['Uniform','High Risk','High Contact'])
redux['Strategy'] = [next(strats) for num in range(12)]
redux['Objective'] = ['Infections']*3 + ['Hospitalizations']*3 + ['Deaths']*3 + ['YLLs']*3
redux['% Reduction'] = list(percent_reduc_infs[0:3])+list(percent_reduc_hosps[0:3])+list(percent_reduc_deaths[0:3])+list(percent_reduc_YLL[0:3])
redux['Number of Individuals'] = list(infs[0:3])+list(hosps[0:3])+list(deaths[0:3])+list(YLLs[0:3])

sns.catplot(x='Objective',y='% Reduction',hue = 'Strategy', data = redux, kind = 'bar')
redux.to_excel('USACaseStudy.xlsx')
#ax = sns.histplot(df, x='Objective', hue='Strategy', weights='% Reduction',multiple='Type', palette='tab20c', shrink=0.8)
#ax.set_ylabel('percentage')
# Fix the legend so it's not on top of the bars.
#legend = ax.get_legend()
#legend.set_bbox_to_anchor((1, 1))
plt.ylim(50, 100)
plt.show()
#sns.catplot(x='Objective',y='Number of Individuals',hue = 'Strategy', data = redux, kind = 'bar',sharey = False)
#plt.show()
rcumul_infs
