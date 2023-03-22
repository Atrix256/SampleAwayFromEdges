import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd

graphs = [
    {
        "title": "Regular Smooth",
        "source": "out/1DResultsSmooth.csv",
        "dest": "out/1DSmoothRegular",
        "cols":[
            "White",
            "Stratified",
            "Golden Ratio",
            "Regular - Ends",
            "Regular - Left",
            "Regular - Center",
            "Regular - Center Equal",
            ]
    },
    {
        "title": "Regular Non Smooth",
        "source": "out/1DResultsNonSmooth.csv",
        "dest": "out/1DNonSmoothRegular",
        "cols":[
            "White",
            "Stratified",
            "Golden Ratio",
            "Regular - Ends",
            "Regular - Left",
            "Regular - Center",
            "Regular - Center Equal",
            ]
    },
    {
        "title": "Blue Smooth",
        "source": "out/1DResultsSmooth.csv",
        "dest": "out/1DSmoothBlue",
        "cols":[
            "White",
            "Stratified",
            "Golden Ratio",
            "Blue - Wrap",
            "Blue - No Wrap",
            "Blue - No Wrap Edge",
            "Blue - No Wrap Half Edge",
            ]
    },        
    {
        "title": "Blue Non Smooth",
        "source": "out/1DResultsNonSmooth.csv",
        "dest": "out/1DNonSmoothBlue",
        "cols":[
            "White",
            "Stratified",
            "Golden Ratio",
            "Blue - Wrap",
            "Blue - No Wrap",
            "Blue - No Wrap Edge",
            "Blue - No Wrap Half Edge",
            ]
    },    
]

for graph in graphs:
    print(graph["title"])
    df = pd.read_csv(graph["source"])

    fig, ax = plt.subplots(1, 1, figsize=(16, 12))

    ax.set_title(graph["title"])

    for col in graph["cols"]:
        line, = ax.plot(df[col])
        line.set_label(col)

    ax.legend()

    ax.set_xscale('log', base=2)
    ax.set_yscale('log', base=2)

    ax.set_xlim(1, 200)

    fig.tight_layout()
    fig.savefig(graph["dest"] + ".png", bbox_inches='tight')
    fig.savefig(graph["dest"] + ".pdf", bbox_inches='tight')

# Setup a plot such that only the bottom spine is shown
def SetupNumberline(ax, title):
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(which='major', width=1.00)
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', width=0.75)
    ax.tick_params(which='minor', length=2.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.patch.set_alpha(0.0)

    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_minor_locator(ticker.NullLocator())

    ax.text(0.0, 0.2, title, transform=ax.transAxes)    

    #ax.xaxis.set_major_locator(ticker.AutoLocator())
    #ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())    

# Make numberlines

cols = [
    "Regular - Ends",
    "Regular - Left",
    "Regular - Center",
    "Regular - Center Equal",
    "Golden Ratio",
    "Blue - Wrap",
    "Blue - No Wrap",
    "Blue - No Wrap Edge",
    "Blue - No Wrap Half Edge",
    "Stratified",
    "White",
]

df = pd.read_csv("out/1DPoints.csv")

fig, ax = plt.subplots(len(cols), 1, figsize=(8, 6))

fig.suptitle("Numberlines")

for colIndex in range(len(cols)):
    SetupNumberline(ax[colIndex], cols[colIndex])
    for x in df[cols[colIndex]]:
        ax[colIndex].plot(x, 0, 'ro', ms = 3, mfc = 'r', clip_on=False, zorder=100)

fig.tight_layout()
fig.savefig("out/1Dnumberline.png", bbox_inches='tight')
fig.savefig("out/1Dnumberline.pdf", bbox_inches='tight')