import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

A = pd.read_csv("3328_WGS.tsv", sep="\t")

Ag = A.groupby("Chromosome")
chr_ends = Ag["End.bp"].apply(max).cumsum()
chr_ends[0] = 0
chr_ends = chr_ends.sort_index()

plt.figure(1); plt.clf()
ax = plt.axes()
for chrom, C in Ag:
    for _, x in C.iterrows():
        span = x[["Start.bp", "End.bp"]] + chr_ends[chrom - 1]

        # plot mus
        ax.plot(
          span,
          x["mu.major"]*np.r_[1, 1],
          color = 'r',
          linestyle = (0, [1, 1]),
          linewidth = 3
        )
        ax.plot(
          span,
          x["mu.minor"]*np.r_[1, 1],
          color = 'b',
          linestyle = (1, [1, 1]),
          linewidth = 3
        )
        # draw solid line across each dotted segment for clarity
        ax.plot(
          span,
          x["mu.major"]*np.r_[1, 1],
          color = 'r',
          linewidth = 0.5
        )
        ax.plot(
          span,
          x["mu.minor"]*np.r_[1, 1],
          color = 'b',
          linewidth = 0.5
        )

        # plot sigmas
        ax.add_patch(patches.Rectangle(
          [span[0], x["mu.major"] - x["sigma.major"]],
          span[1] - span[0], 2*x["sigma.major"],
          color = "r",
          alpha = 0.5,
          linewidth = 0
        ))
        ax.add_patch(patches.Rectangle(
          [span[0], x["mu.minor"] - x["sigma.minor"]],
          span[1] - span[0], 2*x["sigma.minor"],
          color = "b",
          alpha = 0.5,
          linewidth = 0
        ))

    ax.axvline(chr_ends[chrom], linestyle = ":", color = 'k')

ax.set_xlim([0, chr_ends.iloc[-1]])
plt.xticks((chr_ends.iloc[1:].values + chr_ends[:-1].values)/2, chr_ends.iloc[1:].index)