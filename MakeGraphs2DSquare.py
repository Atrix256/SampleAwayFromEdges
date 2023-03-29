import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import math

graphs = [
    {
        "title": "Regular Smooth",
        "source": "out/2DSquareResultsSmooth.csv",
        "dest": "out/2DSquareSmoothRegular",
        "cols":[
            "White",
            "Stratified",
            "Regular Grid",
            "Hex Grid",
            ]
    },
    {
        "title": "Regular Non Smooth",
        "source": "out/2DSquareResultsNonSmooth.csv",
        "dest": "out/2DSquareNonSmoothRegular",
        "cols":[
            "White",
            "Stratified",
            "Regular Grid",
            "Hex Grid",
            ]
    },
    {
        "title": "Regular Non Separable",
        "source": "out/2DSquareResultsNonSeparable.csv",
        "dest": "out/2DSquareNonSeparableRegular",
        "cols":[
            "White",
            "Stratified",
            "Regular Grid",
            "Hex Grid",
            ]
    },
    {
        "title": "LDS Smooth",
        "source": "out/2DSquareResultsSmooth.csv",
        "dest": "out/2DSquareSmoothLDS",
        "cols":[
            "White",
            "Stratified",
            "R2",
            "Halton23",
            "Burley Sobol",
            "Sobol",
            "Fibonacci",
            ]
    },
    {
        "title": "LDS Non Smooth",
        "source": "out/2DSquareResultsNonSmooth.csv",
        "dest": "out/2DSquareNonSmoothLDS",
        "cols":[
            "White",
            "Stratified",
            "R2",
            "Halton23",
            "Burley Sobol",
            "Sobol",
            "Fibonacci",
            ]
    },
    {
        "title": "LDS Non Separable",
        "source": "out/2DSquareResultsNonSeparable.csv",
        "dest": "out/2DSquareNonSeparableLDS",
        "cols":[
            "White",
            "Stratified",
            "R2",
            "Halton23",
            "Burley Sobol",
            "Sobol",
            "Fibonacci",
            ]
    },
    {
        "title": "Blue Smooth",
        "source": "out/2DSquareResultsSmooth.csv",
        "dest": "out/2DSquareSmoothBlue",
        "cols":[
            "White",
            "Stratified",
            "Blue - Wrap",
            "Blue - No Wrap",
            "Blue - No Wrap Edge",
            "Blue - No Wrap Half Edge",
            ]
    },        
    {
        "title": "Blue Non Smooth",
        "source": "out/2DSquareResultsNonSmooth.csv",
        "dest": "out/2DSquareNonSmoothBlue",
        "cols":[
            "White",
            "Stratified",
            "Blue - Wrap",
            "Blue - No Wrap",
            "Blue - No Wrap Edge",
            "Blue - No Wrap Half Edge",
            ]
    },
    {
        "title": "Blue Non Separable",
        "source": "out/2DSquareResultsNonSeparable.csv",
        "dest": "out/2DSquareNonSeparableBlue",
        "cols":[
            "White",
            "Stratified",
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

def SetupPointPlot(ax, title):

    if len(title) == 0:
        ax.spines['right'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')

    #ax.spines['right'].set_color('none')
    #ax.spines['left'].set_color('none')
    #ax.yaxis.set_major_locator(ticker.NullLocator())
    #ax.spines['top'].set_color('none')
    #ax.xaxis.set_ticks_position('bottom')
    #ax.tick_params(which='major', width=1.00)
    #ax.tick_params(which='major', length=5)
    #ax.tick_params(which='minor', width=0.75)
    #ax.tick_params(which='minor', length=2.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.patch.set_alpha(0.0)

    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.xaxis.set_minor_locator(ticker.NullLocator())

    ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_minor_locator(ticker.NullLocator())

    ax.set_title(title)
    #ax.text(0.0, 0.2, title, transform=ax.transAxes)    

    #ax.xaxis.set_major_locator(ticker.AutoLocator())
    #ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())    

# Make numberlines

cols = [
    "White",
    "Stratified",
    "Regular Grid",
    "Hex Grid",
    "R2",
    "Halton23",
    "Sobol",
    "Burley Sobol",
    "Fibonacci",
    "Blue - Wrap",
    "Blue - No Wrap",
    "Blue - No Wrap Edge",
    "Blue - No Wrap Half Edge",
]

numCols = 3
numRows = int(math.ceil(len(cols) / numCols))

df = pd.read_csv("out/2DPoints.csv")

fig, ax = plt.subplots(numRows, numCols, figsize=(4 * numCols, 4 * numRows))

#fig.suptitle("Points")

for colIndex in range(len(cols)):
    axisX = colIndex % numCols
    axisY = int(colIndex / numCols)
    SetupPointPlot(ax[axisY][axisX], cols[colIndex])
    for x,y in zip(df[cols[colIndex] + " X"], df[cols[colIndex] + " Y"]):
        ax[axisY][axisX].plot(x, y, 'ro', ms = 3, mfc = 'r', clip_on=False, zorder=100)

if (len(cols)) % numCols > 0:
    for axisX in range((len(cols)) % numCols, numCols):
        SetupPointPlot(ax[numRows-1][axisX], "")

fig.tight_layout()
fig.savefig("out/2DSquarePoints.png", bbox_inches='tight')
fig.savefig("out/2DSquarePoints.pdf", bbox_inches='tight')