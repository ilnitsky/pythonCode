
import numpy as np
import pandas as pd
from pandas import DataFrame, read_csv
from astropy.table import Table, Column
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

import json
import subprocess
import requests
import pandas as pd
import io

# query = 'piwil1'
query = '63152at7147'
# query = 'Protein prickle'
level = '7147'
species = '7147'

with open('Diptera.csv') as organisms_list:
    organisms_df = read_csv('Diptera.csv')


with open("OG.csv") as OGS:
    OG_list = read_csv('OG.csv', sep=';')['label']

l = []
for query in OG_list:
    # url = 'https://v101.orthodb.org/tab?query=161067at50557&level=50557&species=50557&universal=&singlecopy='
    tab = 'https://v101.orthodb.org/tab?query=%s&level=%s&species=2759&universal=&singlecopy=' % (query, level)
    # fasta = 'https://v101.orthodb.org/fasta?query=%s&level=%s&species=%s&universal=&singlecopy=' % (query, level, species)

    tabData = requests.get(tab).content
    # fastaData = requests.get(fasta).content


    rawData = pd.read_csv(io.StringIO(tabData.decode('utf-8')),sep='\t')
    
    og_name = rawData.groupby(['og_name', 'pub_og_id']).size().reset_index(name='prot_count')['og_name'][0]
    print(og_name)
    rawData = rawData[['organism_name', 'pub_gene_id']]
    for i in organisms_df['organism']:
        if i in rawData['organism_name'].values.tolist():
            l.append(1)
        else:
            l.append(0)
        # print(l)
    organisms_df[og_name] = pd.DataFrame(l,columns=[og_name])
    l = []

organisms_df = organisms_df.set_index('organism')   
# print(organisms_df)

new_df = organisms_df
Z = linkage(new_df.T)
old_order_list = new_df.T.index.tolist()
new_order_list = leaves_list(Z).tolist() 
new_order = [old_order_list[i] for i in new_order_list]
dendrogram(Z)
new_df = new_df[new_order]
new_df = new_df.iloc[::-1]
new_df.to_csv('saved.csv')
print(new_df)


fig = go.Figure(data = go.Heatmap(
    z = new_df.to_numpy(),
    x = new_df.columns.values.tolist(),
    y = new_df.index.tolist(),
    xgap = 2,
    ygap = 2, 
    colorscale = 'Blues'
))
fig.show()


