import pandas as pd

import sys

filename = sys.argv[1]

f = open(filename, 'r')

data = {}
for line in f:
    els = line.split(',')
    print(els)
    label = els[0]
    vals = els[1:-1]
    if label not in data:
        data[label] = []
    data[label] += [float(v) for v in vals]

print("Data:", data)

df = pd.DataFrame(data)

print(df)

means = df.mean()

print(means)
