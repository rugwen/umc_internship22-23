# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 19:17:00 2022

@author: gkoen
"""
#important notes! remove FFPE.csv and change the file directory
import pandas as pd
from scipy.stats import entropy
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from math import log
import matplotlib.patches as mpatches

resp_color = "#1b4ccf" #you can change the color
nonresp_color = "#cf1b1b"

blue_patch = mpatches.Patch(color=resp_color, label='Responders')
red_patch = mpatches.Patch(color=nonresp_color, label='Non-responders')


# Make list of all patient/sample names by scanning through the data




location = r"/data"

files_to_scan = []
for file in os.listdir(location):
    if ".xlsx" in file:
        pass
    else:
        files_to_scan.append(file)
        
names_of_patients = []


for file in files_to_scan:
    directory = location + r"/" + file
    data = pd.read_csv(directory, delimiter = ";")
    names = list(data.columns)[5:]
    
    for patient in names:
        numberpatient = patient.split("r_")[1].split("_")[0]
        print(numberpatient)
        names_of_patients.append(numberpatient)
    
#%%
    
top = 50

for sample in names_of_patients:
    curdf = pd.DataFrame()
    for file in files_to_scan:
        directory = location + r"/" + file
        data = pd.read_csv(directory, delimiter = ";")
        names = list(data.columns)[5:]
        for patient in list(data.columns)[5:]:
            print(patient)

            
            treatment = patient.split("_GK")[0][-1]
            
            
            if sample in patient.split(".")[1]:
                print(sample, patient)
                
                if "_r_" in patient:
                    responsive = True
                else:
                    responsive = False
                
                data.sort_values(by = [patient], ascending = False, inplace=True)
                data.drop(data.index[top:], inplace = True)
                
                if curdf.empty:
                    curdf["Clonotype"] = data['part2']
                    curdf["Sequence"] = data['part3']
                    curdf = curdf.reset_index(drop=True)
                    
                # zoek de frequenties van de sequences op
                for i in range(len(curdf)):
                    sequence = curdf.iloc[[i]]['Sequence']
                    clonotype = curdf.iloc[[i]]['Clonotype']
                    currow = data[data['part3']==sequence.to_string(index=False)]
                    correctrow = currow[currow['part2']==clonotype.to_string(index=False)]
                    
                    value = correctrow[patient]
                    
                    
                    if value.empty:
                        curdf.at[i, "c"+treatment] = 0
                    else:
                        value = sum(value)
                        curdf.at[i, "c"+treatment] = float(value)

            
            else:
                pass

    if responsive:            
        file_name = fr"data/{sample}_r_top{top}.xlsx"
    elif responsive == False:
        file_name = fr"data/{sample}_nr_top{top}.xlsx"
    curdf.to_excel(file_name) 
                

#%% Plotjes making
#remove patient folder in topf 082 and 901/991

normalize = True # zet op False als je niet genormalizeerde plotjes wil

location = r"data/"

resultsdf = pd.DataFrame(columns = ["Patient", "Responsive", "Cure", "Frequency"])

topn = 1 #change based on interest e.g. top 0 or top 1
cur_treatment = "c4"

files_to_scan = []
for file in os.listdir(location):
    if ".xlsx" in file:
        files_to_scan.append(file)
    else:
        pass
    

for file in files_to_scan:
    directory = location + r"/" + file
    data = pd.read_excel(directory)
    print(file)
    data.sort_values(by = cur_treatment, ascending = False, inplace=True)
    print(data.loc[0])
    if "nr" in file:
        responsive = False
    else:
        responsive = True
    if normalize:
        patientname = file.split("_")[0]
        resultsdf.loc[len(resultsdf)] = [patientname, responsive, "1", data.iloc[topn]['c1']/data.iloc[topn]['c1']] #you change this what you want to compare
        resultsdf.loc[len(resultsdf)] = [patientname, responsive, "2", data.iloc[topn]['c2']/data.iloc[topn]['c1']]
        resultsdf.loc[len(resultsdf)] = [patientname, responsive, "4", data.iloc[topn]['c4']/data.iloc[topn]['c1']]
    else:
        patientname = file.split("_")[0]
        resultsdf.iloc[len(resultsdf)] = [patientname, responsive, "1", data.iloc[topn]['c1']] #this can also be changed based on your interest
        resultsdf.iloc[len(resultsdf)] = [patientname, responsive, "2", data.iloc[topn]['c2']]
        resultsdf.iloc[len(resultsdf)] = [patientname, responsive, "4", data.iloc[topn]['c4']]
        
if normalize:
    file_name = rf"data/top{topn}_c1_across_patients_norm_Blanca_change.xlsx"
else:
    file_name = rf"data/top{topn}_c1_across_patients.xlsx"
    
resultsdf['Frequency'] = np.log2(resultsdf['Frequency']).fillna(0)
resultsdf.replace([np.inf, -np.inf], np.nan, inplace=True)
#resultsdf.dropna(inplace=True)
resultsdf.to_excel(file_name) 

sns.lineplot(x = 'Cure', y = 'Frequency', data=resultsdf, hue="Responsive", palette=[resp_color, nonresp_color], hue_order = [True, False])
plt.legend(handles=[blue_patch, red_patch], loc = "best")
plt.ylabel("log2 fold change relative to c1")
if normalize: 
    plt.savefig( fr"data/top{topn}_{cur_treatment}_dominant_clonotype_fig1_norm_log.png")
else:
      plt.ylabel("Frequency %")
      plt.savefig( fr"data/top{topn}_{cur_treatment}_dominant_clonotype_fig1.png")

plt.show()
plt.clf()


x = ['1', '2', '4']
fig = plt.figure()
ax = plt.axes()
label_y_list = []
label_x_list = []

#https://www.python-graph-gallery.com/web-line-chart-with-labels-at-line-end
#https://stackoverflow.com/questions/16992038/inline-labels-in-matplotlib
for patient in list(set(resultsdf["Patient"])):
    if sum(resultsdf[resultsdf["Patient"] == patient]["Responsive"] == True) > 1:
        ccolor = resp_color
        clabel = "Responder"
    else:
        ccolor = nonresp_color
        clabel = "Non-responder"
    

    
    ax.plot(x, resultsdf[resultsdf["Patient"] == patient]["Frequency"], color = ccolor)
    label_y_pos = resultsdf[resultsdf["Patient"] == patient]["Frequency"].tolist()[2]
    print(label_y_pos, patient)
    #label_x_pos = 2.05
    label_x_pos = 2.05
    for i in label_y_list:
        if label_y_pos*0.95 <= i <=label_y_pos*1.06:
            if i>label_y_pos:
                label_y_pos -=0.01
            if i < label_y_pos:
                label_y_pos +=0.1

    label_y_list.append(label_y_pos)
    label_x_list.append(label_x_pos)
    ax.annotate(patient, (label_x_pos, label_y_pos), color = ccolor, size = 7)

    
    plt.legend(handles=[blue_patch, red_patch], loc = "best")
plt.xlim(0, 2.2)
plt.xlabel("Cure")
    
plt.show
if normalize: 
    plt.ylabel(f"log2 fold change relative to c1")
    plt.savefig( fr"data/top{topn}_{cur_treatment}_dominant_clonotype_fig2_norm_change.svg")
else: 
    plt.ylabel("Frequency %")
    plt.savefig( fr"data/top{topn}_{cur_treatment}_dominant_clonotype_fig2.svg")
