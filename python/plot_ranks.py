import numpy as np
import matplotlib.pyplot as plt


# prg seed, rank p_s , rank p, rank q

def process_data(filename):
    data = np.genfromtxt(filename,delimiter=',', dtype=np.float64)
    data = data[:,1:]

    # create difference of rank P - rank Q column
    diff = data[:,1]-data[:,2]
    data = np.insert(data,3,diff,axis=1)
    mean = np.mean(data,axis=0)
    low = np.min(data, axis=0)
    high = np.max(data, axis=0)
   
    
    
    return (mean, low, high)

# NOTE change the director depending on g: g0 or g1.
def accumulate_files(n):
    meanarr = []
    lowarr = [] 
    higharr = [] 
    for m in range(40,301,10):
        # for g=0
        filename = f"g0/data{m}x{n}"
        (mean, low, high) = process_data(filename)
        meanarr += [mean]
        lowarr += [low]
        higharr += [high]
    meanarr = np.array(meanarr)
    lowarr = np.array(lowarr)
    higharr = np.array(higharr)
    mvalues = np.arange(40,301,10) 
    return (meanarr, lowarr, higharr, mvalues)

def plot_rank(mean, low, high, xaxis, labels, loc = "lower right"):
    
    fig, ax = plt.subplots(1)
    ax.plot(xaxis, mean, lw=1, label="mean", color='black', ls='--')
    ax.fill_between(xaxis, low, high, facecolor='blue', alpha=0.3,
                label='range')

    ax.legend(loc=loc)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.grid()

    plt.title(labels[2])
    plt.show()


# for n = 50
(mean, low, high, xaxis) = accumulate_files(50)
# plot rank of P

plot_rank(mean[:,0], low[:,0], high[:,0], xaxis, ["m (number of rows in P)", "rank" , "Rank of $P_\mathbf{s}$ when n = 50"])
#plot_rank(mean[:,1], low[:,1], high[:,1], xaxis, ["m (number of rows in P)", "rank" , "Rank of P when n = 50"])
#plot_rank(mean[:,2], low[:,2], high[:,2], xaxis, ["m (number of rows in P)", "rank" , "Rank of Q when n = 50"])
#plot_rank(mean[:,3], low[:,3], high[:,3], xaxis, ["m (number of rows in P)", "rank P-rank Q" , "Difference between ranks of P and Q when n = 50"], "upper right")
