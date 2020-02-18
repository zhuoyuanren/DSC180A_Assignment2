import pandas as pd
import json
import os

def get_table(team, year):
    df = pd.read_html('https://www.pro-football-reference.com/teams/' + team.lower() + '/' + str(year) + '.htm')
    return df[1]

def get_data(team,year,outpath,**kwargs):
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    for t in team:
        for y in year:
            table = get_table(t, y)
            path = ("%s/%s_%s.csv"%(outpath,t,str(y)))
            table.to_csv(path)
