import pandas as pd
import os

pd.set_option("display.line_width", 300)

df = pd.read("../../data/cfa/cfasnIa_param.dat", index_col="sndata")

dfx.describe()