import numpy as np
import matplotlib.pyplot as plt

# prg seed, status (0 success, 1 ran out of attempts, -1 P not full rank, 2 wrong string found), #attempts

def process_data(filename, off):
    data = np.genfromtxt(filename,delimiter=',', dtype=np.float64)
    #strip the prg seed
    data = off-data[:,1]

    mean = np.mean(data)
    std = np.std(data)
    low = mean-std
    high = mean+std
   
    return (mean, low, high)

def accumulate_files(dire, mul):
    meanarr = []
    lowarr = [] 
    higharr = []

    for n in range(40,301,10):
        # for g=0
        filename = f"{dire}/data{n}"
        (mean, low, high) = process_data(filename, min(n,n*mul))

        # ignore tests where everything was skipped
        meanarr += [mean]
        lowarr += [low]
        higharr += [high]

    meanarr = np.array(meanarr)
    lowarr = np.array(lowarr)
    higharr = np.array(higharr)
    nvalues = np.arange(40,301,10) 

    return (meanarr, lowarr, higharr, nvalues)


def plot_mean(ax, mean, low, high, xaxis):

    ax.plot(xaxis, mean, lw=1, label="mean", color='black', ls='--')
    ax.fill_between(xaxis, low, high, facecolor='blue', alpha=0.3,
                label='range within 1 stddev')

  #  ax.legend(loc="upper left")
    ax.set_xlabel("n")
    ax.set_ylabel("min(n,m) - rank $M^TM$")
    ax.grid()


fig, ax = plt.subplots(nrows=1, ncols=3)


(mean, low, high, xaxis) = accumulate_files("rand_P",2)
plot_mean(ax[0], mean, low, high, xaxis)
ax[0].set_title("m = 2n")
(mean, low, high, xaxis) = accumulate_files("rand_P2",0.5)
plot_mean(ax[1], mean, low, high, xaxis)
ax[1].set_title("m=n/2")
(mean, low, high, xaxis) = accumulate_files("rand_P3",1)
plot_mean(ax[2], mean, low, high, xaxis)
ax[2].set_title("m=n")

handles, labels = ax[2].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right')

fig.suptitle("Rank of $M^TM$ for a random M")
fig.tight_layout()
plt.show()
