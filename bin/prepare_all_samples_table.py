#!/usr/bin/env python3

import json
import sys
import pandas as pd

jsons = sys.argv[1]

data = []
for path in jsons.split(" "):
    with open(path) as f:
        data.append(json.load(f))

df = pd.DataFrame(data)
df.to_csv("all_samples.csv", index=False)
