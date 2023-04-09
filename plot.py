from plot_tree import draw_tree
from plot_msa import my_ax
import matplotlib.pyplot as plt

# Start with a square Figure.
fig = plt.figure(figsize=(12, 3))
# Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
# the size of the marginal axes and the main axes in both directions.
# Also adjust the subplot parameters for a square plot.
gs = fig.add_gridspec(1, 2,  width_ratios=(1, 3),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)
# Create the Axes.
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

# Throught the get_ylim() function to acquire the ymin and ymax

my_ax("example.fas",ax2,"example.nwk",ymin=0.2,ymax=6.8).get_ax()
draw_tree("example.nwk",ax=ax1)



