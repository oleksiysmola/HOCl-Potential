import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)
# from decimal import Decimal, getcontext

# getcontext().prec = 20


# columns = ["grep", "rCH1", "rCH2", "rCH3", "Sa", "Sb", "rho", "E", "point"]
# columns = ["grep", "rCH", "rCO", "alpha", "E", "E_VQZ", "point"]
# columns = ["grep", "rCH", "rCO", "alpha", "E", "E2", "E3", "point"]
columns = ["grep", "rCH", "rCO", "alpha", "E", "E2", "point"]
# df = pd.read_csv("CH3OH-3DEnergies.dat", delim_whitespace=True, names=columns, dtype=str)
df = pd.read_csv("HOCl-MRCI-MoreAngle.en", delim_whitespace=True, names=columns, dtype=str)
# df2 = pd.read_csv("HOCl-VQZ-Grid.en", delim_whitespace=True, names=["point", "E_VQZ"], dtype=str)
df = df[["rCH", "rCO", "alpha", "E", "point"]]
# df2 = pd.read_csv("VQZ.energies", delim_whitespace=True, names=["point", "E_VQZ"], dtype=str)
print(df.head(10).to_string(index=False))
df["point"] = df["point"].astype(float).astype(int)
# df2["point"] = df2["point"].astype(float).astype(int)
# def splitPoint(row):
#     row["point"] = row["point"].split(".")[0]
#     return row
print(df.head(10).to_string(index=False))
# df = df.parallel_apply(lambda x: splitPoint(x), axis=1, result_type="expand")
# df2 = df2.parallel_apply(lambda x: splitPoint(x), axis=1, result_type="expand")
# df = pd.merge(df, df2, "left", on="point")
df["point"] = df["point"].astype(int)
# df["E_value"] = df["E_VQZ"].astype(float)
df["E_value"]  = df["E"].astype(float)*2.194746354e05
df["E_value"] = df["E_value"] - df["E_value"].min()
df = df.sort_values(by="point")
print(df[df["E_value"] >= 50000].to_string(index=False))
df = df[df["rCH"].astype(float) <= 2.3]
df = df[df["rCO"].astype(float) <= 2.25]
df = df[df["E_value"] <= 50000]
# df = pd.merge(df, df2, on="point", how="left")
# df = df[["rCH", "rCO", "alpha", "E", "E_VQZ", "point"]]
# df2 = df[["rCH", "rCO", "alpha", "E_VQZ", "point"]]
df = df.to_string(index=False, header=False)
# df2 = df2.to_string(index=False, header=False)
# print(df)
vqzFile = "HOCl_MRCI-MoreAngle.dat"
# cbsFile = "HOCl_CBS.dat"
with open(vqzFile, "w+") as FileToWriteTo:
    FileToWriteTo.write(df)
# with open(cbsFile, "w+") as FileToWriteTo:
#     FileToWriteTo.write(df1)
# re1= 1.4296
# re2= 0.95887
# re3= 1.092294
# re4= 1.092294
# re5= 1.092294
# ae1= 107.9812
# ae2= 110.6646
# ae3= 110.6646
# ae4= 110.6646
# rco   =re1+ s1
# rch0  =re2+ s2
# rch1  =re3+ s3
# rch2  =re4+ s4
# rch3  =re5+ s5
# acoh  =ae1+ s6
# aoch1 =ae2+ s7
# aoch2 =ae2+ s8
# aoch3 =ae2+ s9