import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import math

graphs = [
    {
        "title": "Regular Smooth",
        "source": "out/2DCircleResultsSmooth.csv",
        "dest": "out/2DCircleSmoothRegular",
        "cols":[
            "White",
            "Stratified",
            "Fibonacci",            
            "Regular Grid",
            "Hex Grid",
            ]
    },
    {
        "title": "Regular Non Smooth",
        "source": "out/2DCircleResultsNonSmooth.csv",
        "dest": "out/2DCircleNonSmoothRegular",
        "cols":[
            "White",
            "Stratified",
            "Fibonacci",            
            "Regular Grid",
            "Hex Grid",
            ]
    },
    {
        "title": "Regular Non Separable",
        "source": "out/2DCircleResultsNonSeparable.csv",
        "dest": "out/2DCircleNonSeparableRegular",
        "cols":[
            "White",
            "Stratified",
            "Fibonacci",
            "Regular Grid",
            "Hex Grid",
            ]
    },    
    {
        "title": "LDS Smooth",
        "source": "out/2DCircleResultsSmooth.csv",
        "dest": "out/2DCircleSmoothLDS",
        "cols":[
            "White",
            "Stratified",
            "Fibonacci",
            "R2",
            "Halton23",
            "Halton23 Circle",
            "Sobol",
            "Sobol Circle",
            "Burley Sobol",
            "Burley Sobol Circle",
            ]
    },
    {
        "title": "LDS Non Smooth",
        "source": "out/2DCircleResultsNonSmooth.csv",
        "dest": "out/2DCircleNonSmoothLDS",
        "cols":[
            "White",
            "Stratified",
            "Fibonacci",            
            "R2",
            "Halton23",
            "Halton23 Circle",
            "Sobol",
            "Sobol Circle",
            "Burley Sobol",
            "Burley Sobol Circle",
            ]
    },
    {
        "title": "LDS Non Separable",
        "source": "out/2DCircleResultsNonSeparable.csv",
        "dest": "out/2DCircleNonSeparableLDS",
        "cols":[
            "White",
            "Stratified",
            "Fibonacci",            
            "R2",
            "Halton23",
            "Halton23 Circle",
            "Sobol",
            "Sobol Circle",
            "Burley Sobol",
            "Burley Sobol Circle",
            ]
    },    
    {
        "title": "Blue Smooth",
        "source": "out/2DCircleResultsSmooth.csv",
        "dest": "out/2DCircleSmoothBlue",
        "cols":[
            "White",
            "Stratified",
            "Blue - No Wrap",
            "Blue - No Wrap Edge",
            "Blue - No Wrap Half Edge",
            ]
    },        
    {
        "title": "Blue Non Smooth",
        "source": "out/2DCircleResultsNonSmooth.csv",
        "dest": "out/2DCircleNonSmoothBlue",
        "cols":[
            "White",
            "Stratified",
            "Blue - No Wrap",
            "Blue - No Wrap Edge",
            "Blue - No Wrap Half Edge",
            ]
    },
    {
        "title": "Blue Non Separable",
        "source": "out/2DCircleResultsNonSeparable.csv",
        "dest": "out/2DCircleNonSeparableBlue",
        "cols":[
            "White",
            "Stratified",
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
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    #ax.yaxis.set_major_locator(ticker.NullLocator())
    ax.spines['top'].set_color('none')
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

    if len(title) > 0:
        circle = plt.Circle((0.5, 0.5), 0.5, color='b', fill=False, clip_on=False)
        ax.add_patch(circle)

    ax.set_title(title)
    #ax.text(0.0, 0.2, title, transform=ax.transAxes)    

    #ax.xaxis.set_major_locator(ticker.AutoLocator())
    #ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())    

# Make numberlines

cols = [
    "White",
    "Regular Grid",
    "Regular Grid Circle",
    "Stratified",
    "Stratified Circle",
    "Hex Grid",
    "Hex Grid Circle",
    "R2",
    "R2 Circle",
    "Halton23",
    "Halton23 Circle",
    "Sobol",
    "Sobol Circle",
    "Burley Sobol",
    "Burley Sobol Circle",
    "Fibonacci",
    "Blue - No Wrap",
    "Blue - No Wrap Edge",
    "Blue - No Wrap Half Edge",
]

numCols = 4
numRows = int(math.ceil(len(cols) / numCols))

df = pd.read_csv("out/2DCirclePoints.csv")

fig, ax = plt.subplots(numRows, numCols, figsize=(4 * numCols, 4 * numRows), squeeze=False)

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
fig.savefig("out/2DCirclePoints.png", bbox_inches='tight')
fig.savefig("out/2DCirclePoints.pdf", bbox_inches='tight')

# Make the secondary numberline

cols = [
    "Fibonacci Simple",
    "Regular Lattice Circle"
]

fig, ax = plt.subplots(1, 3, figsize=(12, 4))
SetupPointPlot(ax[0], "Fibonacci (Red) vs Lattice (Green)")
SetupPointPlot(ax[1], "Fibonacci")
SetupPointPlot(ax[2], "Lattice")

for x,y in zip(df[cols[0] + " X"], df[cols[0] + " Y"]):
    ax[0].plot(x, y, 'ro', ms = 3, mfc = 'r', clip_on=False, zorder=100)
    ax[1].plot(x, y, 'ro', ms = 3, mfc = 'r', clip_on=False, zorder=100)

for x,y in zip(df[cols[1] + " X"], df[cols[1] + " Y"]):
    ax[0].plot(x, y, 'go', ms = 3, mfc = 'g', clip_on=False, zorder=100)
    ax[2].plot(x, y, 'go', ms = 3, mfc = 'g', clip_on=False, zorder=100)

fig.tight_layout()
fig.savefig("out/2DCirclePointsB.png", bbox_inches='tight')
fig.savefig("out/2DCirclePointsB.pdf", bbox_inches='tight')
