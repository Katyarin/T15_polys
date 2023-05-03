import json
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import requests
from pathlib import Path
import txt_to_json

shots = [42090]
polys = [34, 35]
G = 10
las_energy = 1.6

plot_osc = False
plot_exp_sig = False # False
time_for_start_plot = 1 #ms

print(shots)
#polyn = 34
thomson = 'ust' #or 'usual' #or 'ust' or 'divertor'

sp_cal = '22.06.30'
abs_cal = '22.06.01'
#laser_const = 625448756.1404866 #881940473.5373197
#laser_const = 485448756.1404866 #881940473.5373197

bound1_prof = 160
bound2_prof = 180
ne_corr_mult_all = 1
ne_corr_mult_T15 = 13.2e20 #8.33e40
delta_time = 3.73

URL = 'https://172.16.12.130:443/api'
TSpath = 'TS_core/'
raw_data_path = 'data/DRS/plasma/raw_data/'
res_data_path = 'data/DRS/plasma/results/'

with open('source/%s_abs_cal.json' %abs_cal, 'r') as abs_file:
    abs_data = json.load(abs_file)

A = abs_data['A']


#const
sigma_ts = 8/3 * math.pi * 2.8179403267 * 1e-15 * 2.8179403267 * 1e-15
lamd0 = 1064 *1e-9
q = 1.6e-19 #cl

const_for_ne = A * lamd0 * sigma_ts / q
print(const_for_ne)


def TS_res(dev, N_phe, sigma, i, err):
    ne = []
    chi = []
    for j in range(len(dev['ch']['1'])):
        num = 0
        den = 0
        ch_count = 0
        for ch in dev['ch'].keys():
            if err[ch][i] == 'off scale' or N_phe[ch][i] < 0:
                continue
            num += N_phe[ch][i] * dev['ch'][ch][j] / math.pow(sigma[ch][i], 2)
            den += dev['ch'][ch][j] * dev['ch'][ch][j] / math.pow(sigma[ch][i], 2)
            ch_count+=1
        if ch_count < 2:
            return 1, 1, 1000, 100, 100
        ne_loc = num / den
        ne.append(num / den)
        chi_local = 0
        for ch in dev['ch'].keys():
            if N_phe[ch][i] == None:
                continue
            chi_local += (N_phe[ch][i] - ne_loc * dev['ch'][ch][j]) ** 2 / math.pow(sigma[ch][i],2)
        chi.append(chi_local)
    chi_min = min(chi)
    index = chi.index(chi_min)
    if index == len(dev['Te']) - 1:
        index = len(dev['Te']) - 2
        #print('index out of range')
    Te = dev['Te'][chi.index(chi_min)]
    ne_res = ne[chi.index(chi_min)]

    '''sigma_Te'''
    dN_dTe = {}
    for ch in dev['ch'].keys():
        dN_dTe[ch] = np.diff([i * ne_res for i in dev['ch'][ch]]) / np.diff(dev['Te'])
        #dN_dTe[ch].append(dev['ch'][ch][-1] * ne_res)
    sum1 = 0
    sum2 = 0
    sum3 = 0
    for ch in dev['ch'].keys():
        if N_phe[ch][i] == None:
            continue
        sum1 += (dev['ch'][ch][index] * ne_res / sigma[ch][i]) ** 2
        sum2 += (dN_dTe[ch][index] / sigma[ch][i]) ** 2
        sum3 += dN_dTe[ch][index] / sigma[ch][i] * dev['ch'][ch][index] * ne_res / sigma[ch][i]
    sum3 = sum3 ** 2
    sigma_Te = (sum1 / (sum1 * sum2 - sum3)) ** 0.5
    wtf: float = 1e-40
    sigma_ne: float = 1e40
    if (sum1 * sum2 - sum3) != 0 and (las_energy * const_for_ne) != 0:
        sigma_ne = ne_res / (las_energy * const_for_ne) * (sum2 / (sum1 * sum2 - sum3)) ** 0.5
        wtf = ne_res / (las_energy * const_for_ne)
    #(i, ne_res, N_phe['6'][i])

    return Te, wtf, chi_min, sigma_Te, sigma_ne


def find_start_integration(signal):
    maximum = signal.index(max(signal))
    for i in range(0, maximum):
        if signal[maximum - i] > 0 and signal[maximum - i - 1] <= 0:
            return maximum - i - 1
    return 0


def find_end_integration(signal):
    maximum = signal.index(max(signal))
    for i in range(maximum, len(signal) - 1):
        if signal[i] > 0 and signal[i + 1] <= 0:
            return i + 1
    return len(signal) - 1


def to_phe(shotn):

    filename = res_data_path + '%d/%d.json' %(shotn, shotn)

    M = 100
    el_charge = 1.6 * 10 ** (-19)
    #G = 10
    R_sv = 10000
    freq = 5  # GS/s
    time_step = 1 / freq  # nanoseconds
    event_len = 1024

    delta = {0: 0, 1: 150, 2: 170, 3: 190, 4: 200, 5: 210, 6: 0}
    timestamps = {}
    N_photo_el = {}
    var_phe = {}
    timeline = [i * time_step for i in range(event_len)]
    calc_err = {}
    #print(timeline)
    for polyn in polys:
        N_photo_el[polyn] = {}
        var_phe[polyn] = {}
        calc_err[polyn] = {}
        timestamps[polyn] = []

    for polyn in polys:
        for ch in range(7):
            N_photo_el[polyn][ch] = []
            var_phe[polyn][ch] = []
            calc_err[polyn][ch] = []


    with open(filename, 'r') as file:
        raw_data0 = json.load(file)
    for polyn in polys:
        raw_data = raw_data0[str(polyn)]
        p = 0
        for event in raw_data:
            timestamps[polyn].append(event['t']/1000 + delta_time)

            #if max(event['ch'][1]) > 0.030:
            if event['t']/1000 + delta_time > time_for_start_plot and plot_osc:
                fig2, axs2 = plt.subplots(3, 3)
                fig2.suptitle(event['t']/1000 + delta_time)
                p = 1
            for ch in range(7):
                signal = event['ch'][ch]
                if max(signal) > 0.8:
                    calc_err[polyn][ch].append('off scale')
                    #print(event['t']/1000 + delta_time)
                else:
                    calc_err[polyn][ch].append('')
                if ch == 0:
                    pre_sig = 100
                else:
                    pre_sig = 200
                base_line = sum(signal[0:pre_sig]) / len(signal[0:pre_sig])
                for i in range(len(signal)):
                    signal[i] = signal[i] - base_line

                if ch == 0:
                    index_0 = 0
                    for i, s in enumerate(signal[10:]):
                        if s > 0.250:
                            index_0 = i - 20
                            #print(index_0)
                            break

                for i in range(len(signal)):
                    signal[i] = signal[i] * 1000
                var_in_sr = np.var(signal[0:pre_sig])
                delta_exp = {}
                if thomson == 'usual':
                    width = 120
                    delta_exp[ch] = delta[ch]
                elif thomson == 'ust':
                    width = 270
                    delta_exp[ch] = delta[ch] - 100
                elif thomson == 'divertor':
                    width = 100
                    delta_exp[ch] = delta_divertor[ch]
                else:
                    print('something wrong! Unnown config')
                    stop
                start_index = index_0 + delta_exp[ch]
                end_index = start_index + width


                if p:
                    #print(max(event['ch'][5]), event['t']/1000)
                    axs2[int(ch//3), int(ch%3)].set_title('ch = ' + str(ch))
                    axs2[int(ch//3), int(ch%3)].plot(signal)
                    axs2[int(ch//3), int(ch%3)].vlines(start_index, min(signal), max(signal))
                    axs2[int(ch//3), int(ch%3)].vlines(end_index, min(signal), max(signal))


                Ni = np.trapz(signal[start_index:end_index],
                                    timeline[start_index:end_index]) / (M * el_charge * G * R_sv * 0.5)
                if ch==6:
                    Ni = Ni * (M * el_charge * G * R_sv * 0.5)
                N_photo_el[polyn][ch].append(Ni *1e-12)
                var = math.sqrt(math.fabs(6715 * 0.0625 * var_in_sr - 1.14e4 * 0.0625) + math.fabs(Ni *1e-12) * 4)
                var_phe[polyn][ch].append(var)
            p=0
        plt.figure(figsize=(10, 3))
        plt.title('Shot #' + str(shotn) + ', poly #' + str(polyn))
        for ch in N_photo_el[polyn].keys():
            #color = ['r', 'g', 'b', 'm', 'black', 'orange', 'brown', 'pink']
            if ch != 0 and ch != 6:
                plt.errorbar(timestamps[polyn], N_photo_el[polyn][ch], yerr=var_phe[polyn][ch], label='ch' + str(ch))
                plt.scatter([t for i, t in enumerate(timestamps[polyn]) if calc_err[polyn][ch][i] == 'off scale'],
                            [j for i, j in enumerate(N_photo_el[polyn][ch]) if calc_err[polyn][ch][i] == 'off scale'], marker='x', s=40, c='black', zorder=2.5)
            #plt.plot(timestamps[polyn], N_photo_el[polyn][ch], '^-', label='ch' + str(ch))
        N_photo_el[polyn][6][0] = N_photo_el[polyn][6][1]
        plt.ylabel('N, phe')
        plt.grid()
        plt.xlabel('time')
        plt.legend()
    for_temp = {}
    for polyn in polys:
        for_temp[polyn] = {'timeline': timestamps[polyn], 'data': N_photo_el[polyn], 'err': var_phe[polyn], 'culc_err': calc_err[polyn],
                    'laser_en': N_photo_el[polyn][6]}
    with open(res_data_path + '%d/N_phe.json' %shotn, 'w') as f:
        json.dump(for_temp, f)


def Temp(shotn):

    '''data in phe'''
    with open(res_data_path + '%d/N_phe.json' %shotn, 'r') as fp:
        data = json.load(fp)

    for_temp = {}

    for polyn in polys:
        '''waiting signals'''
        with open('source/' + sp_cal + '_dev_num_' + str(polyn) + '.json', 'r') as file:
            dev_num = json.load(file)
        polyn = str(polyn)
        N_phe = data[polyn]['data']
        timeline = data[polyn]['timeline']
        sigma = data[polyn]['err']
        culc_errors = data[polyn]['culc_err']

        Te = []
        Te_err = []
        ne = []
        ne_err = []
        chi2 = []
        start_temp = 18
        end_temp = 65
        for i in range(len(N_phe['0'])):
            # for i in range(start_temp, end_temp):
            Te_loc, ne_loc, chi2_loc, Te_err_loc, ne_err_loc = TS_res(dev_num, N_phe, sigma, i, culc_errors)
            if plot_exp_sig and timeline[i] > time_for_start_plot:
            #if Te_err_loc/Te_loc < 0.2:
                plt.figure()
                plt.title(timeline[i])
                for ch in dev_num['ch'].keys():
                    plt.plot(dev_num['Te'], [i * ne_loc * (las_energy * const_for_ne) for i in dev_num['ch'][ch]], label=ch)
                    if culc_errors[ch][i] == 'off scale':
                        plt.scatter(Te_loc, N_phe[ch][i], marker='x')
                    else:
                        plt.scatter(Te_loc, N_phe[ch][i])
                #plt.ylim(0, max([i * ne_loc for i in dev_num['ch']['4']]) * 1.25)
                plt.legend()
                plt.show()
            Te.append(Te_loc)
            ne.append(ne_loc)
            chi2.append(chi2_loc)
            Te_err.append(Te_err_loc)
            ne_err.append(ne_err_loc)
        for_temp[polyn] = {'timeline': timeline, 'Te': Te, 'Te_err': Te_err, 'ne': ne, 'ne_err': ne_err, 'chi2': chi2 }


    with open(res_data_path + '%d/Te.json' %shotn, 'w') as res_file:
        json.dump(for_temp, res_file)

    return for_temp


"""___________________________________________________________________________________________________________________"""

T15 = {'Te': [], 'Te_err': [], 'ne': [], 'ne_err': [], 'chi2': []}
p = 0
for shotn in shots:
    txt_to_json.to_json(shotn, save_file=True)

    to_phe(shotn)

    figTe, axTe = plt.subplots(1,1)
    figne, axne = plt.subplots(1, 1)
    for polyn in polys:
        T15_data = Temp(shotn)[str(polyn)]
        print(T15_data.keys())
        print(len(T15_data['timeline']))
        #print(len(data['data']['events']))

        axTe.errorbar([t for i, t in enumerate(T15_data['timeline']) if T15_data['Te_err'][i] / T15_data['Te'][i] < 0.7],
                     [t for i, t in enumerate(T15_data['Te']) if T15_data['Te_err'][i] / T15_data['Te'][i] < 0.7],
                     yerr=[t for i, t in enumerate(T15_data['Te_err']) if T15_data['Te_err'][i] / T15_data['Te'][i] < 0.7], ls='--', label=polyn)
        axTe.set_ylim(0, 2000)
        plt.grid()
        axTe.set_xlabel('time, ms')
        axTe.set_ylabel('Te, eV')
        axTe.legend()

        axne.errorbar([t for i, t in enumerate(T15_data['timeline']) if T15_data['ne_err'][i] / T15_data['ne'][i] < 0.7],
                     [t*ne_corr_mult_T15 for i, t in enumerate(T15_data['ne']) if T15_data['ne_err'][i] / T15_data['ne'][i] < 0.7],
                     yerr=[t*ne_corr_mult_T15 for i, t in enumerate(T15_data['ne_err']) if T15_data['ne_err'][i] / T15_data['ne'][i] < 0.7], ls='--', label=polyn)
        #plt.ylim(0, 2000)
        axne.grid()
        axne.set_xlabel('time, ms')
        axne.set_ylabel('ne')
        axne.legend()
        #plt.show()

        with open(res_data_path + '%d/%i_result_dynamic.txt' %(shotn, polyn), 'w') as f_res1:
            f_res1.write(' %14s' % 'time_ms')
            f_res1.write(' %14s' % 'Te')
            f_res1.write(' %14s' % 'Te_err')
            f_res1.write(' %14s' % 'ne')
            f_res1.write(' %14s' % 'ne_err')
            f_res1.write(' %14s' % 'chi2')
            f_res1.write('\n')
            for i in range(len(T15_data['timeline'])):
                f_res1.write(' %14.4f' % T15_data['timeline'][i])
                f_res1.write(' %14.4f' % T15_data['Te'][i])
                f_res1.write(' %14.4f' % T15_data['Te_err'][i])
                f_res1.write(' %14.4f' % (T15_data['ne'][i]*ne_corr_mult_T15))
                f_res1.write(' %14.4f' % (T15_data['ne_err'][i]*ne_corr_mult_T15))
                f_res1.write(' %14.4f' % T15_data['chi2'][i])
                f_res1.write('\n')

plt.show()
