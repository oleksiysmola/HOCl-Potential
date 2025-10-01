import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

columns = ["point", "DBOC"]

df = pd.read_csv("greppedDBOC-CL35_aug-cc-pVTZ+d.txt", delim_whitespace=True, names=columns, dtype=str)

# from decimal import Decimal, getcontext

# getcontext().prec = 20
df = df.drop_duplicates()

def splitPoint(row):
    row["point"] = row["point"].split(".")[0].split("_")[-1]
    return row
df = df.parallel_apply(lambda x: splitPoint(x), axis=1, result_type="expand")
# df = df[["point", "rOH", "rOCL", "A", "E"]]
# df = df[["rCO", "rOH", "rCH1", "rCH2", "rCH3", "aHOC","aOCH1", "aOCH2", "raOCH3", "dH1", "dH2", "dH3", "E", "point"]]
# df[["dH1", "dH2", "dH3"]] = (df[["dH1", "dH2", "dH3"]].astype(float) + 120.00)
df["point"] = df["point"].astype(int)
# displacements = df[["dH1", "dH2", "dH3"]]
# displacements["dH1"] = displacements["dH1"] - 180.0
# displacements["dH2"] = displacements["dH2"] - 300.0
# displacements["dH3"] = displacements["dH3"] - 420.0
# df["dH1"] = 180.0 + displacements["dH2"]
# df["dH2"] = 300.0 + displacements["dH3"]
# df["dH3"] = 60.0 + displacements["dH1"]
df = df.sort_values(by="point")

grid = pd.read_csv("/scratch/vol1/asmola/HOCl/molpro/HOCL-GRID.txt", delim_whitespace=True, names=["point", "rOH", "rOCL", "A", "E"], dtype=str)
grid["point"] = grid["point"].astype(int)
df = pd.merge(grid, df, on="point", how="left")
energies = pd.read_csv("CFOUR-Cl35.energies", delim_whitespace=True, names=["point", "Ecfour"], dtype=str)
energies["point"] = energies["point"].astype(int)
df = pd.merge(df, energies, on="point", how="left")
# df = df.to_string(index=False, header=False, formatters={'dH1': '{:.8f}'.format, 'dH2': '{:.8f}'.format, 'dH3': '{:.8f}'.format})
df = df.to_string(index=False, header=False)
statesFile = "HOCL35-SortedDBOCs_aug-cc-pVTZ+d.txt"
with open(statesFile, "w+") as FileToWriteTo:
    FileToWriteTo.write(df)