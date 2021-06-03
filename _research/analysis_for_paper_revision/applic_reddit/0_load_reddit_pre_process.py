#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\applic_reddit\0_load_reddit_pre_process.py
#Created Date: Tuesday May 4th 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Thursday May 6th 2021 1:46:42 pm
#Modified By:  Domenico Di Gangi
#-----
#Description: preprocess reddit hyperlink data downloaded from https://snap.stanford.edu/data/soc-RedditHyperlinks.html
#-----
########



#%%
import pandas as pd
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

#%% 
# load data and rename columns
data_path = "../../../data/reddit_hyperlinks/raw_data/"
os.listdir(data_path)
col_names = ["source", "target", "post_id", "time", "sentiment", "properties"]  
df_orig = pd.read_csv(f"{data_path}soc-redditHyperlinks-body.tsv", names = col_names,  sep="\t", header = 0)

df_orig["datetime"] = pd.to_datetime(df_orig.time)
df_orig = df_orig.set_index("datetime")
df_orig = df_orig.sort_values(by="datetime")

#%% EDA

# check aggregate number of obs
df_count = df_orig.time.resample("W").count()
plt.plot(df_count, ".")
plt.plot(df_count[df_count==0], ".r")
# number of nodes appearing at least once
pd.concat((df_orig.source, df_orig.target)).unique().shape[0]
# how many nodes appear n-times ? 
plt.plot(pd.concat((df_orig.source, df_orig.target)).value_counts().value_counts(), ".")
plt.yscale("log")
plt.xscale("log")
#%% EDA
# check daily statistics over time

df = df_orig

df_aggr_stats = df.resample("W").agg(lambda x :pd.Series({
    "n_nodes" : len(np.unique(x[["source", "target"]].values.ravel("K"))), 
    "n_links" : len(x), 
    "out_l_95q" : x["source"].value_counts().quantile(0.95),
    "in_l_95q" : x["target"].value_counts().quantile(0.95),
    "out_l_99q" : x["source"].value_counts().quantile(0.99),
    "in_l_99q" : x["target"].value_counts().quantile(0.99),
     } ) )


df_aggr_stats.plot(y= ["out_l_95q", "in_l_95q", "out_l_99q", "in_l_99q"], title="quantiles of degree distributions")
plt.figure()
(df_aggr_stats.n_links / (df_aggr_stats.n_nodes * (df_aggr_stats.n_nodes-1)) ).plot(title="density", logy=True)



# %%
# save only columns required for discrete time modeling
df = df_orig
df = df[["source", "target"]].reset_index()
df["week"] = df["datetime"].dt.isocalendar().week
df["date"] = df["datetime"].dt.date
df["hour"] = df["datetime"].dt.hour

df[["date", "hour", "source", "target"]].to_csv(f"{data_path}reddit_hyperlinks_edges.csv")

# %%
