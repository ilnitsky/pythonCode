
import numpy as np
import pandas as pd
from pandas import DataFrame, read_csv
from astropy.table import Table, Column
import json
import subprocess
import requests
import pandas as pd
import io
import time

start_time = time.time()

level = '2759'
species = '9606'


with open("OG.csv") as OGS:
    OG_list = read_csv('OG.csv', sep=';')['label']

j = 0
l = []
for query in OG_list:
    json_relevant = 'https://www.orthodb.org/pgrest/rpc/search?query=%s&level=%s&skip=0&limit=100' % (query, level)
    # print(json_relevant)
    jsonData = requests.get(json_relevant).content
    jsonData = json.loads(jsonData)
    print(jsonData["data"], jsonData["bigdata"][0]["name"])
    og_json = jsonData["data"][0]

    tab = 'https://v101.orthodb.org/tab?query=%s&level=%s&species=%s&universal=&singlecopy=' % (og_json, level, species)    
    tabData = requests.get(tab).content

    rawData = pd.read_csv(io.StringIO(tabData.decode('utf-8')),sep='\t')
    
    print(rawData["pub_gene_id"].tolist())  
    j += 1
    print(str(j)+"/"+str(len(OG_list)))



print("--- %s seconds ---" % (time.time() - start_time))
