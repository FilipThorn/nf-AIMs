#!/usr/bin/env python 

#Modified script by Filip Thörn. Original from jupyter-notebook written by Mozes Blom

try:
	import os
	import sys
	import subprocess
	import pandas as pd
	import datetime
	import argparse

except ImportError:
	sys.exit("One of the required modules can't be found...")

parser = argparse.ArgumentParser()
parser.add_argument("-F", "--Parents_fst", help="WEIR_AND_COCKERHAM_FST for parental populations", action = "store", required=True)
parser.add_argument("-I", "--P1H_fst", help="WEIR_AND_COCKERHAM_FST for parent1 and Hybrid", action = "store", required=True)
parser.add_argument("-j", "--P2H_fst", help="WEIR_AND_COCKERHAM_FST for parent2 and Hybrid", action = "store", required=True)
parser.add_argument("-P", "--P1", help="Parent 1", action = "store", required=True)
parser.add_argument("-p", "--P2", help="Parent 2", action = "store", required=True)
parser.add_argument("-O", "--Out", help="Output location and name", action = "store", required=True)
args = parser.parse_args()


#varibles
Parents_fst = args.Parents_fst
Hyb_parent1_fst = args.P1H_fst
Hyb_parent2_fst = args.P2H_fst
Parent1 = args.P1
Parent2 = args.P2
Output_table = args.Out



print('AIMs calculation started: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
sys.stdout.flush()

####read data####

#read parental Fst file
parents_fst_df = pd.read_csv(Parents_fst,
                             header = 0, sep="\t")

# Only retain sites with Fst == 1
parents_fixed_df = parents_fst_df[parents_fst_df['WEIR_AND_COCKERHAM_FST'] == 1]

# Import the Fst estimates between the hybrid parent1
Hyb_parent1_df = pd.read_csv(Hyb_parent1_fst,
                             header = 0, sep="\t")
# Import the Fst estimates between the hybrid parent2
Hyb_parent2_df = pd.read_csv(Hyb_parent2_fst,
                             header = 0, sep="\t")

print('Import data done: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

#Intersect the fixed sites, with the Fst sites for the hybrid:
## Please note this weird concat call was used to avoid the set copy warning error!
parents_fixed_df = pd.concat([parents_fixed_df, (parents_fixed_df["CHROM"] + '_' + parents_fixed_df["POS"].astype(str)).rename("combo")], axis=1)
Hyb_parent1_df = pd.concat([Hyb_parent1_df, (Hyb_parent1_df["CHROM"] + '_' + Hyb_parent1_df["POS"].astype(str)).rename("combo")], axis=1)
Hyb_parent2_df = pd.concat([Hyb_parent2_df, (Hyb_parent2_df["CHROM"] + '_' + Hyb_parent2_df["POS"].astype(str)).rename("combo")], axis=1)

merged_parent1_df = parents_fixed_df.merge(Hyb_parent1_df, left_on='combo', right_on='combo', indicator=True)
merged_parent2_df = parents_fixed_df.merge(Hyb_parent2_df, left_on='combo', right_on='combo', indicator=True)

print('Joining done: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

# Remove excess columns and rename specific columns
merged_parent1_df = merged_parent1_df.drop(["CHROM_x", "POS_x", "WEIR_AND_COCKERHAM_FST_x", "CHROM_y", "POS_y", "_merge"], axis=1)
merged_parent2_df = merged_parent2_df.drop(["CHROM_x", "POS_x", "WEIR_AND_COCKERHAM_FST_x", "CHROM_y", "POS_y", "_merge"], axis=1)
merged_parent1_df = merged_parent1_df.rename(columns={"WEIR_AND_COCKERHAM_FST_y": Parent1})
merged_parent2_df = merged_parent2_df.rename(columns={"WEIR_AND_COCKERHAM_FST_y": Parent2})


# merge the two processed dataframes, should be equal in size so no need for indicator
hybrid_df = merged_parent1_df.merge(merged_parent2_df, left_on='combo', right_on='combo')


# Replace NaN with 0
hybrid_df.fillna(0, inplace = True)


# Function that infers the genotype relative to the parental species (see above text)
def f(row):
    if (row[Parent1] != 0) and (row[Parent2] != 0):
        if (abs(row[Parent1] - row[Parent2]) < 0.2):
            val = Parent1 + 'X' + Parent2
        else:
            val = 'NaN'
    elif (row[Parent1] == 0) and (row[Parent2] == 1):
        val = Parent1 + 'X' + Parent1
    elif (row[Parent1] == 1) and (row[Parent2] == 0):
        val = Parent2 + 'X' + Parent2
    else:
        val = 'NaN'
    return val

hybrid_df["geno"] = hybrid_df.apply(f, axis=1)

# remove the NaN individuals and split the combo column
# Match the chromosome and position information (for downstream plotting)
hybrid_df = hybrid_df[hybrid_df.geno != 'NaN']

hybrid_df[['chromo','position']] = hybrid_df.combo.str.split("_",expand=True)

print('Begin printing to csv: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))

hybrid_df.to_csv(Output_table, index = False, header=True)

print('AIMs calculation finished: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
sys.stdout.flush()
