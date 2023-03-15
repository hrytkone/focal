import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap
import sys
import ROOT

np.set_printoptions(precision=1)

if len(sys.argv) != 2:
    print("USAGE : %s <input file>" % (sys.argv[0]))
    sys.exit(1)

infilename = sys.argv[1]
fin = ROOT.TFile.Open(infilename, "READ")

# Get histograms
histPt = fin.Get("hRecPtVsTruePt")
histEta = fin.Get("hRecEtaVsTrueEta")
histE = fin.Get("hRecEVsTrueE")

fig = plt.figure()
gs = fig.add_gridspec(3, hspace=0)
axs = gs.subplots(sharex=True, sharey=True)
fig.suptitle('Sharing both axes')
axs[0].hist2d(histPt, cmap=plt.cm.jet)
axs[1].plot(histEta)
axs[2].plot(histE)
plt.show()

# Hide x labels and tick labels for all but bottom plot.
#for ax in axs:
#    ax.label_outer()