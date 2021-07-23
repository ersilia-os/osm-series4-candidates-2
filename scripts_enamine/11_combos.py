from __init__ import OUTPUT

import numpy as np
import pandas as pd
import os

df = pd.read_csv(os.path.join(OUTPUT, "data_10.csv"))

print(df.shape)

df.to_csv(os.path.join(OUTPUT, "data_11.csv"), index = False)
