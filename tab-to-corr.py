
from operator import not_
import numpy as np
import pandas as pd
from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

import pandas as pd
import time



#Load the XML Species and Taxids

with open('vertebrata-taxid.xml', 'r') as file:
    lines = file.readlines()
    i = 0
    Names = []
    Taxids = []
    while i < (len(lines)-1):
        if "<id>"  in lines[i+1]:
            Names.append(lines[i].lstrip(' ').rstrip(' ').lstrip('<name>').rstrip('\n').rstrip('</name>'))
            Taxids.append(lines[i+1].lstrip(' ').rstrip(' ').lstrip('<id>').rstrip('\n').rstrip('</id>'))

        i += 1

    df_txs = pd.DataFrame({
        'Names': Names,
        'Taxids': Taxids
        }) 
     
    Taxid_list = df_txs['Taxids'].tolist()

level = "Vertebrata_7742.tab"
df = pd.read_csv(level, delimiter="\t", header=None)
df.columns = ['OG', 'level', 'species', 'gene']

df[['species', 'm']] = df['species'].str.split('_', expand=True)
print(df)
series = df.groupby(['OG', 'species']).size()

#Filter unused indexes
not_index = ['unknown_0', '302047', '40157', '8962', '8845', '8801', '9708', '9694', '9767', '9901', '9807', 'unknown_1', '9778', '9818', 'unknown_2', 'unknown_3', 'unknown_4', 'unknown_5', 'unknown_6', 'unknown_7', 'unknown_8', 'unknown_9', 'unknown_10', 'unknown_11', 'unknown_12', 'unknown_13', 'unknown_14', 'unknown_15', 'unknown_16', 'unknown_17', 'unknown_18', 'unknown_19', 'unknown_20', 'unknown_21', 'unknown_22', 'unknown_23', 'unknown_24', 'unknown_25', 'unknown_26', 'unknown_27', 'unknown_28']
Taxid_list = [x for x in Taxid_list if x not in not_index]

new_df = series.unstack(level=1).fillna(0)

new_df[new_df != 0] = 1
new_df = new_df[Taxid_list]




# new_df = new_df.iloc[0:1000,:]
new_df = new_df.T.add_suffix('at7742').T
new_df['91919191919'] = 0.0



start_time = time.time()


#Compute Pearson Correlation

corr = new_df.T.corr()

# display(corr)
corr.to_csv("corr.csv")


#Compute top correlations
diff_l = {}
diff_s = {}

for column in corr.columns:
    diff_l[column] = corr[column].nlargest(11)[1:].index.tolist()
    diff_s[column] = corr[column].nsmallest(10).index.tolist()
    
print(len(diff_l))

diff_lx = pd.DataFrame(diff_l).T
diff_sx = pd.DataFrame(diff_s).T
diff_slx = pd.concat([diff_lx, diff_sx], axis=1)

# print(diff_slx)
diff_slx.to_csv("diff_lsx.csv")

print("--- %s seconds ---" % (time.time() - start_time))



