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
marvelEnergies = pd.read_csv("HOCl-MARVEL.en", delim_whitespace=True, names=columns)
marvelEnergies = marvelEnergies.sort_values(by="J")

marvelEnergies = marvelEnergies.to_string(index=False)

marvelFile = "HOCl-Marvel-JSort.en"
with open(marvelFile, "w+") as FileToWriteTo:
    FileToWriteTo.write(marvelEnergies)
