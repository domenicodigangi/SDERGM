#########
#File: c:\Users\digan\Dropbox\Dynamic_Networks\repos\ScoreDrivenExponentialRandomGraphs\_research\analysis_for_paper_revision\applic_collegMsg\0_load_college_msg_pre_process.py
#Created Date: Friday April 23rd 2021
#Author: Domenico Di Gangi,  <digangidomenico@gmail.com>
#-----
#Last Modified: Tuesday May 4th 2021 11:51:20 am
#Modified By:  Domenico Di Gangi
#-----
#Description: Preprocess college messages data downloaded from SNAP :https://snap.stanford.edu/data/CollegeMsg.html
#-----
########





#%%
import pandas as pd
import numpy as np
import os
import sys
from matplotlib import pyplot as plt

#%% 
# load data
data_path = "../../../data/collegeMsg/raw_data/"
os.listdir(data_path)
col_names = ["source", "target", "time"] 
{n: pd.UInt32Dtype for n in col_names}
df_orig = pd.read_csv(f"{data_path}CollegeMsg.txt", names = col_names, dtype  ={n: np.uint64() for n in col_names}, sep=" ")

df_orig["datetime"] = pd.to_datetime(df_orig.time, unit="s")
df_orig = df_orig.set_index("datetime")
df_orig = df_orig.sort_values(by="datetime")

#%% 
# check daily number of obs
df_count = df_orig.time.resample("D").count()
plt.plot(df_count, ".")
plt.plot(df_count[df_count==0], ".r")

df = df_orig
#%% 
# check daily statistics over time


df_aggr_stats = df.resample("0.5D").agg(lambda x :pd.Series({
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

df = df[["source", "target"]].reset_index()
df["date"] = df["datetime"].dt.date
df["hour"] = df["datetime"].dt.hour

df[["date", "hour", "source", "target"]].to_csv(f"{data_path}collegeMsg_daily_edges.csv")

# %%
