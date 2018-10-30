import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

# TODO ==== add clean up data
# TODO ==== add generation of plot title from data e.g launcher, engine, engine number, departure year, parameters to plot
# TODO ==== automatically generate title
# TODO ==== fix colorbar ticks so they are central

def plot_pareto_power(outbound = True,
                 file = r'C:\Users\andrew mcsweeney\Documents\projects\merlot\examples\merlot_outbound_data.csv',
                      title='merlot_outbound_data'):
    # read data to pandas data frame
    data = pd.read_csv(file)
    df = pd.DataFrame(data, columns=['helio_tof', 'Power', 'orb_arr_mass'])

    # data to plot
    x = df['helio_tof'].values
    y = df['orb_arr_mass'].values
    c = df['Power'].values


    # create colourmap with number of unique points
    # cmap = plt.get_cmap('jet',(max(c)-min(c))+1)
    cmap = plt.get_cmap('jet', len(set(c))-1)
    cmap = plt.get_cmap('jet', 10)

    # create plot
    fig , ax = plt.subplots()
    cax = ax.scatter(x=x, y=y, c=c, s=10, cmap=cmap)
    fig.colorbar(cax)

    # axis limits
    plt.xlim([x.min()-50,x.max()+50])
    # plt.ylim(0,round(y.max() / 1000.0) * 1000.0)

    # axis labels
    plt.xlabel('Total time of flight (days)')
    plt.ylabel('Delivered mass (kg)')

    # plot title
    plt.suptitle(title)

    # show plot
    plt.show()



plot_pareto_power(file = r'C:\Users\andrew mcsweeney\Documents\projects\merlot\run_data\outbound-Mars-2026-ARM-1-falcon9.csv',
                  title = 'outbound-Mars-2026-ARM-1-falcon9')

