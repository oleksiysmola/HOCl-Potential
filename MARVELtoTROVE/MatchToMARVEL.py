import pandas as pd
import numpy as np
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)
from decimal import Decimal, getcontext

getcontext().prec = 20


# columns = ["grep", "rCH1", "rCH2", "rCH3", "Sa", "Sb", "rho", "E", "point"]
# columns = ["grep", "rCH", "rCO", "alpha", "E", "E_VQZ", "point"]
# columns = ["grep", "rCH", "rCO", "alpha", "E", "E2", "E3", "point"]
columns = ["nu1", "nu2", "nu3", "J", "Ka", "Kc",  "E", "unc", "transitions"]
# df = pd.read_csv("CH3OH-3DEnergies.dat", delim_whitespace=True, names=columns, dtype=str)
marvelEnergies = pd.read_csv("HOCl-MARVEL.en", delim_whitespace=True, names=columns, dtype=str)

def findSymmetry(row):
    if int(row["Kc"]) % 2 == 0:
        row["Gamma"] = "A'"
    else:
        row["Gamma"] = "A\""
    return row

marvelEnergies = marvelEnergies.parallel_apply(lambda x:findSymmetry(x), axis=1, result_type="expand")

troveColumns = ["N", "N2", "J", "Sym", "Obs.", "Calc.", "Obs-Calc.", "Weight", "(", "Sym2", ";", "K", ")", "(2", "Sym3", ";2", "v1", "v2", "v3", ")3", "k2", "V1", "V2", "V3"]
troveEnergies = pd.read_csv("HOCL-Trove.en", delim_whitespace=True, names=troveColumns, dtype=str)
print("\n")
print(len(troveEnergies))
# troveEnergies = troveEnergies[troveEnergies["N"].duplicated() == False]
troveEnergies = troveEnergies.head(int(len(troveEnergies)/101))
print("\n")
print(troveEnergies.head().to_string(index=False))
print("\n")
print(len(troveEnergies))

def findMatchingLevel(row, states):
    matchingStates = states[states["J"] == row["J"]]
    matchingStates = matchingStates[matchingStates["Sym"] == row["Gamma"]]
    matchingStates["Obs-Calc"] = abs(float(row["E"]) - states["Calc."].astype(float))
    matchingStates = matchingStates.sort_values(by="Obs-Calc")
    matchingState = matchingStates.head(1).squeeze()
    row["Nb"] = matchingState["N"]
    row["Etrove"] = matchingState["Calc."]
    row["Obs-Calc"] = matchingState["Obs-Calc"]
    # row["Gamma"] = matchingState["Gamma"]
    row["K"] = matchingState["K"]
    # row["n1"] = matchingState["n1"]
    # row["n2"] = matchingState["n2"]
    # row["n3"] = matchingState["n3"]
    return row

marvelEnergies = marvelEnergies.parallel_apply(lambda x:findMatchingLevel(x, troveEnergies), result_type="expand", axis=1)
# print(marvelEnergies.head(20).to_string(index=False))
marvelEnergies= marvelEnergies.dropna()
refinementEnergies = marvelEnergies.copy()
marvelEnergies = marvelEnergies[["nu1", "nu2", "nu3", "J", "Ka", "Kc",  "E", "unc", "transitions", "Gamma", "Nb", "Etrove", "Obs-Calc"]]
marvelEnergies = marvelEnergies.sort_values(by="Obs-Calc")
marvelEnergies = marvelEnergies.to_string(index=False)
marvelisedFile = "HOCl-Matched.out"
with open(marvelisedFile, "w+") as FileToWriteTo:
    FileToWriteTo.write(marvelEnergies)
print("New states file ready!")


# refinementEnergies = refinementEnergies[refinementEnergies["J"] <= 10]
refinementEnergies["J"] = refinementEnergies["J"].astype(int)
refinementEnergies = refinementEnergies.sort_values(by="J")
refinementEnergies["weight"] = np.tanh(1 - refinementEnergies["unc"].astype(float))*refinementEnergies["transitions"].astype(int)
refinementEnergies = refinementEnergies[["J", "Gamma", "Nb", "E", "Ka", "nu1", "nu2", "nu3", "weight"]]
refinementEnergies = refinementEnergies.to_string(header=False, index=False)
refinementFile = "HOCl-RefinementEnergies.ref"
with open(refinementFile, "w+") as FileToWriteTo:
    FileToWriteTo.write(refinementEnergies)

# marvelEnergiesReduced = marvelEnergiesReduced.to_string(index=False, header=False)
# marvelisedFileReduced = "CS2-Matched-KRot5-REFIT2-AMES-Reduced-NoB.out"
# with open(marvelisedFileReduced, "w+") as FileToWriteTo:
#     FileToWriteTo.write(marvelEnergiesReduced)
# print("New states file ready!")

# print(troveEnergies.to_string(index=False, header=False))
# print(marvelEnergies.to_string(index=False, header=False)) 

# df2 = pd.read_csv("HOCl-VQZ-Grid.en", delim_whitespace=True, names=["point", "E_VQZ"], dtype=str)
# df = df[["rCH", "rCO", "alpha", "E", "E2",  "point"]]
# df2 = pd.read_csv("VQZ.energies", delim_whitespace=True, names=["point", "E_VQZ"], dtype=str)
# print(df.head(10).to_string(index=False))
# df["point"] = df["point"].astype(float).astype(int)
# df2["point"] = df2["point"].astype(float).astype(int)
# def splitPoint(row):
#     row["point"] = row["point"].split(".")[0]
#     return row
# print(df.head(10).to_string(index=False))
# df = df.parallel_apply(lambda x: splitPoint(x), axis=1, result_type="expand")
# df2 = df2.parallel_apply(lambda x: splitPoint(x), axis=1, result_type="expand")
# df = pd.merge(df, df2, "left", on="point")
# df["point"] = df["point"].astype(int)
# df["E_value"] = df["E_VQZ"].astype(float)
# df["E_value"]  = df["E"].astype(float)*2.194746354e05
# df["E_value"] = df["E_value"] - df["E_value"].min()
# df = df.sort_values(by="point")
# print(df[df["E_value"] >= 50000].to_string(index=False))
# df = df[df["rCH"].astype(float) <= 2.3]
# df = df[df["rCO"].astype(float) <= 2.25]
# df = df[df["E_value"] <= 50000]
# df = pd.merge(df, df2, on="point", how="left")
# df = df[["rCH", "rCO", "alpha", "E", "E_VQZ", "point"]]
# df2 = df[["rCH", "rCO", "alpha", "E_VQZ", "point"]]
# df = df.to_string(index=False, header=False)
# df2 = df2.to_string(index=False, header=False)
# print(df)
# vqzFile = "HOCl_MRCI3.dat"
# cbsFile = "HOCl_CBS.dat"
# with open(vqzFile, "w+") as FileToWriteTo:
#     FileToWriteTo.write(df)