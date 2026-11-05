"""One colour scheme for every figure (Okabe-Ito, colour-blind safe). Each model keeps the same colour in every
figure; marker shapes / line styles are kept too, so figures still read correctly when printed in greyscale."""
import matplotlib.pyplot as plt

OURS = "#0072B2"      # per-OH model (this work)       blue
IDEAL = "#7A7A7A"     # ideal (Raoult) baseline         grey
UNIFAC = "#E69F00"    # UNIFAC                          orange
OV = "#009E73"        # osmotic virial                  green
KM = "#CC79A7"        # Kleinhans-Mazur                 purple
DATA = "#1A1A1A"      # measured data                   near-black
ALT = "#D55E00"       # second data set / warning       vermillion
SKY = "#56B4E9"       # light blue
CLASS = {"acyclic polyol": OURS, "sugar": OV, "mono-alcohol": UNIFAC}


def setup(size=7.5):
    plt.rcParams.update({"font.size": size, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "axes.edgecolor": "#333333", "axes.titleweight": "bold", "legend.frameon": False})
