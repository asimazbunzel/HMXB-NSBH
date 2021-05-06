# Matplotlib style sheet
# For more info, see: https://matplotlib.org/users/customizing.html

#### LINES
lines.linewidth      : 1.5     ## line width in points
lines.marker         : None    ## the default marker
lines.markersize     : 4       ## markersize, in points

#### FONT
font.family          : STIXGeneral
font.style           : normal
font.size            : 10.0
font.serif           : DejaVu Serif
font.sans-serif      : DejaVu Sans
font.monospace       : DejaVu Sans Mono

#### TEXT
text.usetex          : False
mathtext.rm          : sans
mathtext.fontset     : stix

#### AXES
axes.labelsize       : 10.0   ## fontsize of the x any y labels

#### TICKS
## see http://matplotlib.org/api/axis_api.html#matplotlib.axis.Tick
xtick.top            : True   ## draw ticks on the top side
xtick.bottom         : True   ## draw ticks on the bottom side
xtick.labelbottom    : True   ## draw label on the bottom
xtick.labelsize      : 10.0   ## fontsize of the tick labels
xtick.direction      : inout  ## direction: in, out, or inout
xtick.minor.visible  : True   ## visibility of minor ticks on x-axis
ytick.left           : True   ## draw ticks on the left side
ytick.right          : True   ## draw ticks on the right side
ytick.labelleft      : True   ## draw tick labels on the left side
ytick.labelsize      : 10.0   ## fontsize of the tick labels
ytick.direction      : inout  ## direction: in, out, or inout
ytick.minor.visible  : True   ## visibility of minor ticks on y-axis

#### Legend
legend.frameon       : False     ## if True, draw the legend on a background patch
legend.fontsize      : 8.0
legend.handlelength  : 2.0      ## the length of the legend lines
legend.handleheight  : 0.7      ## the height of the legend handle
legend.borderpad     : 0.4      ## border whitespace
legend.labelspacing  : 0.2      ## the vertical space between the legend entries
legend.handletextpad : 0.8      ## the space between the legend line and legend text
legend.borderaxespad : 0.7      ## the border between the axes and legend edge

#### FIGURE
## See http://matplotlib.org/api/figure_api.html#matplotlib.figure.Figure
figure.figsize       : 7,6   ## figure size in inches
figure.dpi           : 200

#### SAVING FIGURES
## the default savefig params can be different from the display params
## e.g., you may want a higher resolution, or to make the figure
## background white
savefig.dpi         : figure
savefig.format      : pdf
