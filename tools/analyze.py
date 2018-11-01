
import numpy as np
import matplotlib.pyplot as plt
import pandas

Ndrop = 20

# Main figure
fig, ax = plt.subplots(1, 1)
 
# Filter output
df = pandas.read_csv('../filtered.txt', index_col='t')
df.drop(df.drop(df.index[:Ndrop], inplace=True))
base_t = df.index[0]
df.index -= base_t
df.plot(ax=ax)

# Mocap
dfMocap = pandas.read_csv('../inMocap.txt', index_col='t')
dfMocap.drop(dfMocap.index[:Ndrop], inplace=True)
dfMocap.index -= base_t
dfMocap.plot(ax=ax, style='.')

# IMU
dfIMU = pandas.read_csv('../inIMU.txt', index_col='t')
dfIMU.drop(dfIMU.index[:Ndrop], inplace=True)
dfIMU.index -= base_t
dfIMU.plot()

plt.show()
