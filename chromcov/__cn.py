from pathlib import Path
import pandas as pd

summ = pd.read_csv("data/testmosslow.mosdepth.summary.txt", sep="\t")
autosomes = summ[summ['chrom'].str.match(r"chr\d+\b")]
baseline = autosomes["bases"].sum() / autosomes["length"].sum()

def get_flag(ratio):
	if ratio > 1.25:
		return "gain"
	elif ratio < .75:
		return "loss"
	else:
		return ""

summ = summ[summ['chrom'].str.match(r"chr[\dXY]+")]
summ["ratio"] = summ["mean"] / baseline
summ["approx_cn"] = 2*summ["ratio"]
summ["flag"] = summ["ratio"].apply(lambda r: get_flag(r))

summ.to_csv("out/cn_approx.csv", index=False)