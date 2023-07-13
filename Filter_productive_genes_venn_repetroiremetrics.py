"""
Created on Tue Jul 26 18:36:26 2022

@author: gkoen
"""
#note before using this script. Change file directory.
import pandas as pd
from scipy.stats import entropy
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from math import log
import matplotlib.patches as mpatches
import statannotations
from statannotations.Annotator import Annotator


resp_color = "#1b4ccf" #you can fill in any color you like
nonresp_color = "#cf1b1b"

blue_patch = mpatches.Patch(color=resp_color, label='Responders')
red_patch = mpatches.Patch(color=nonresp_color, label='Non-responders')



location = r"(location)" # r"add file directory"

files_to_scan = []
for file in os.listdir(location):
    if ".xlsx" in file:
        pass  
    else:
        files_to_scan.append(file)
        
        
def chao1_bias_corrected(counts_table):
    """Calculates bias-corrected chao1 given counts: Eq. 2 in EstimateS manual.
    Formula: chao1 = S_obs + N_1(N_1-1)/(2*(N_2+1)) where N_1 and N_2 are
    count of singletons and doubletons respectively.
    Note: this is the bias-corrected formulat from Chao 1987, Eq. 2 in the
    EstimateS manual.
    https://github.com/pycogent/pycogent/blob/master/cogent/maths/stats/alpha_diversity.py
    """
    total = (counts_table!=0).sum()
    singles = (counts_table==1).sum()
    doubles = (counts_table==2).sum()
    return total + singles*(singles-1) / (2.0*(doubles+1))

new_file = location.split(".csv")[0] + "_results.csv"

entropies = []
clonalities = []

resultsdf = pd.DataFrame(columns= ["File", "Study number", "Treatment", "Responsive", "Shannon Entropy", "Pielou Index", "Clonality", "Richness"])
venndf = pd.DataFrame(columns = ["080k1p2", "080k1p3","080k2p2", "080k2p3", "080k3p2","080k3p3","080k4p2","080k4p3","080FFPE", "901-991k1p2", "901-991k1p3","901-991k2p2", "901-991k2p3", "901-991k3p2","901-991k3p3","901-991k4p2","901-991k4p3","901-991FFPE"])

overlapping = ["080", "901-991", "022"]
venndf = pd.DataFrame()

for file in files_to_scan:
    directory = location + r"/" + file
    data = pd.read_csv(directory, delimiter = ";")
    data = data[data['part4'].str.contains('pro') & ~data['part4'].str.contains('unp|pop')]
    names = list(data.columns)[5:]

    for patient in list(data.columns)[5:]: 
        if "_r_" in patient:
            responsive = True
        else:
            responsive = False
            
        if "PaxGene" in patient:
            numberpatient = patient.split("r_")[1].split("_")[0]
            kuur = patient.split("_GK")[0][-1]
            print("Number: ", numberpatient, "Kuur: ", kuur)
            
            
            shannon_entropy = entropy(data[patient][data[patient]!=0].value_counts())
            entropies.append(shannon_entropy)
            print("Shannon entropy:", shannon_entropy)
            
            Hmax = log(len(data[patient][data[patient]!=0]))
            pielou = shannon_entropy/Hmax
            print("Pielou Index:", pielou)
            clonality = 1 - pielou
            
            richness_counts = data[patient][data[patient]!=0]*1/data[patient][data[patient]!=0].min()
            richness_counts = richness_counts.astype(np.int64)
            richness = chao1_bias_corrected(richness_counts)
            
            print("Richness: ", richness, "\n")
            newrow = [file, numberpatient, kuur, responsive, shannon_entropy, pielou, clonality, richness]
            resultsdf.loc[len(resultsdf)] = newrow
            
            
            # with open(new_file, 'w') as f:
            #     writer = csv.writer(f)
            #     writer.writerow(names)
            #     writer.writerow(entropies)
            
            # data collection for Venn diagram patient 80 & 901-991
            if numberpatient in overlapping:
                print("Saving Venn diagram data..")
                colname = numberpatient+"_k_"+kuur
                p2 = data['part2'][data[patient]!=0]
                p3 = data['part3'][data[patient]!=0]
                
                venndf[colname+"_p2"] = p2
                venndf[colname+"_p3"] = p3
                venndf.reset_index(inplace=True,drop=True)
        else:
            print("Skipped", patient)

#FFP uitlezen
for file in files_to_scan:
    directory = location + r"/" + file
    data = pd.read_csv(directory, delimiter = ";")
    data = data[data['part4'].str.contains('pro') & ~data['part4'].str.contains('unp|pop')]
    names = list(data.columns)[5:]

    for patient in list(data.columns)[5:]: 
        if "_r_" in patient:
            responsive = True
        else:
            responsive = False
            
        if "FFPE" in patient:
            numberpatient = patient.split("r_")[1].split("_")[0]
            print("Number: ", numberpatient)
            
            
            shannon_entropy = entropy(data[patient][data[patient]!=0].value_counts())
            entropies.append(shannon_entropy)
            print("Shannon entropy:", shannon_entropy)
            
            Hmax = log(len(data[patient][data[patient]!=0]))
            pielou = shannon_entropy/Hmax
            print("Pielou Index:", pielou)
            clonality = 1 - pielou
            
            richness_counts = data[patient][data[patient]!=0]*1/data[patient][data[patient]!=0].min()
            richness_counts = richness_counts.astype(np.int64)
            richness = chao1_bias_corrected(richness_counts)
            
            print("Richness: ", richness, "\n")
            newrow = [file, numberpatient, kuur, responsive, shannon_entropy, pielou, clonality, richness]
            resultsdf.loc[len(resultsdf)] = newrow
            
            # data collection for Venn diagram patient 80 & 901-991
            if numberpatient in overlapping:
                print(numberpatient)
                print("Saving Venn diagram data..")
                colname = numberpatient+"_FFPE"
                p2 = data['part2'][data[patient]!=0]
                p3 = data['part3'][data[patient]!=0]
                
                venndf[colname+"_p2"] = p2
                venndf[colname+"_p3"] = p3
                venndf.reset_index(inplace=True,drop=True)
                venndf = venndf.apply(lambda x: pd.Series(x.dropna().values)).fillna(' ')

# Specify the name of the excel file
file_name = r"(location)+\NGSdata2105.xlsx" # r"add file directory"

venn_file_name = r"(location)+\venn_data.xlsx" # r"add file directory"

# saving the excelsheet
resultsdf.to_excel(file_name)  
venndf.to_excel(venn_file_name, index=False)  

                
resultsdf[resultsdf["Responsive"]==True] # de dataframe waar responsive == True
resultsdf[resultsdf["Responsive"]==True]["Clonality"] # alleen de clonality waar responsive == true

resultsdf[resultsdf["Treatment"]=='1'] # allen K1

resultsdf[resultsdf["Treatment"]=='1'][resultsdf["Responsive"]==True] # Alleen K1 EN responsive
resultsdf[resultsdf["Treatment"]=='1'][resultsdf["Responsive"]==True]["Clonality"] # Alleen K1 EN responsive

#Richness Chao1



#%%


#pairs=[[("1", True), ("4", True)],
       #[("1", True), ("1", False)],
       #[("1", True), ("2", True)],
       #[("1", True), ("2", False)],
       #[("2", True), ("2", False)],
       #[("4", True), ("4", False)],
       #[("2", True), ("4", True)],
       #[("2", False), ("4", False)]]

pairs = [[("2", True), ("4", True)],
         [("2", False), ("4", False)]]

#Plotjes maken
plt.clf()


sig = True


sns.set_style(style="white")

ax = sns.violinplot(x="Treatment", y="Clonality", hue="Responsive",
                    data=resultsdf, palette = [resp_color, nonresp_color],hue_order= [True, False], dpi=1400)
plt.legend(handles=[blue_patch, red_patch], loc = "best")
fig = ax.get_figure()

if sig:
    annotator = Annotator(ax, pairs, data=resultsdf, x="Treatment", y="Clonality", hue= "Responsive", hue_order= [True, False])
    annotator.configure(test='Mann-Whitney', text_format='star', loc='outside')
    annotator.apply_and_annotate()
    #fig.savefig(location+"Clonality_filter.png",bbox_inches = "tight") # r"add file directory"
plt.show()

plt.clf()


ax2 = sns.violinplot(x="Treatment", y="Shannon Entropy", hue="Responsive",
                    data=resultsdf, palette = [resp_color, nonresp_color],hue_order= [True, False], dpi=1400)
plt.legend(handles=[blue_patch, red_patch], loc = "best")
fig = ax2.get_figure()
if sig:
    annotator = Annotator(ax2, pairs, data=resultsdf, x="Treatment", y="Shannon Entropy", hue= "Responsive", hue_order= [True, False])
    annotator.configure(test='Mann-Whitney', text_format='star', loc='outside')
    annotator.apply_and_annotate()
    #fig.savefig(location+"Shannon_filter.png",bbox_inches = "tight") # r"add file directory"
plt.show()

plt.clf()

#pairs=[[]]

ax3 = sns.violinplot(x="Treatment", y="Richness", hue="Responsive",
                    data=resultsdf, palette = [resp_color, nonresp_color],hue_order= [True, False], dpi=1400)
plt.legend(handles=[blue_patch, red_patch], loc = "best")
fig = ax3.get_figure()

if sig:
    annotator = Annotator(ax3, pairs, data=resultsdf, x="Treatment", y="Richness", hue= "Responsive", hue_order= [True, False])
    annotator.configure(test='Mann-Whitney', text_format='star', loc='outside')
    annotator.apply_and_annotate()
    #fig.savefig(location+"Chao1_filter.png", bbox_inches = "tight") # r"add file directory"
plt.show()
plt.clf()




    
    
#%%
#Venn plotjes maken
import matplotlib
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib.cm import plasma
# https://towardsdatascience.com/how-to-create-and-customize-venn-diagrams-in-python-263555527305
matplotlib.rcParams.update({'font.size': 8.5})
plt.rcParams['savefig.dpi'] = 300

values = np.linspace(0, 1, 4)
colors = [plasma(x) for x in values]
colors1 = [colors[0], colors[1], colors[2]]
colors2 = [colors[3], colors[1], colors[2]]
colors3 = [colors[0], colors[3], colors[1]]
colors4 = [colors[0], colors[3], colors[2]]


for patient in overlapping:
    x = [string for string in list(venndf.columns) if patient in string and 'p2' in string]
    print(x)

    try:
        fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 6))

        # c1 c4 FFPE
        total = len(set(venndf[x[0]]).union(set(venndf[x[2]])).union(set(venndf[x[3]])))
        venn3([set(venndf[x[0]]), set(venndf[x[2]]), set(venndf[x[3]])], set_labels=("c1", "c4", "FFPE"), set_colors=colors1,
              subset_label_formatter=lambda x: str(x), #+ "\n" + f"{x/total:1.0%}" + "",
              ax=axes[0])

        # c2 c4 FFPE
        total = len(set(venndf[x[1]]).union(set(venndf[x[2]])).union(set(venndf[x[3]])))
        venn3([set(venndf[x[1]]), set(venndf[x[2]]), set(venndf[x[3]])], set_labels=("c2", "c4", "FFPE"), set_colors=colors2,
              subset_label_formatter=lambda x: str(x), #+ "\n" + f"{x/total:1.0%}" + "",
              ax=axes[1])

        # c1 c2 c4
        total = len(set(venndf[x[0]]).union(set(venndf[x[1]])).union(set(venndf[x[2]])))
        venn3([set(venndf[x[0]]), set(venndf[x[1]]), set(venndf[x[2]])], set_labels=("c1", "c2", "c4"), set_colors=colors3,
              subset_label_formatter=lambda x: str(x), # + "\n" + f"{x/total:1.0%}" + "",
              ax=axes[2])
        
    except:
        # 901 (c1 c2 FFPE)
        total = len(set(venndf[x[0]]).union(set(venndf[x[1]])).union(set(venndf[x[2]])))
        venn3([set(venndf[x[0]]), set(venndf[x[1]]), set(venndf[x[2]])], set_labels=("c1", "c2", "FFPE"), set_colors=colors4,
              subset_label_formatter=lambda x: str(x), # + "\n" + f"{x/total:1.0%}" + "",
              ax=axes[1])
        
        fig.suptitle(f"Responder {patient}", fontsize=14, y=0.9)
    else:
        fig.suptitle(f"Non-responder {patient}", fontsize=14, y=0.9)

    plt.savefig(fr"(location)+ {patient}_clonotype.png", dpi=300) # r"add file directory"
    plt.savefig(fr"(location)+_{patient}_clonotype.svg", dpi=300) # r"add file directory"

    plt.show()
