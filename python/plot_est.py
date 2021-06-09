import numpy as np
import matplotlib.pyplot as plt


# prg seed, status (0 success, 1 ran out of attempts, -1 P not full rank, 2 wrong string found), #attempts

def process_data(filename):
    data = np.genfromtxt(filename,delimiter=',', dtype=np.float64)
    #strip the prg seed
    data = data[:,1:]

    #count how many with each status
    mean = np.mean(data[:,1],axis=0)
    low = np.min(data[:,1], axis=0)
    high = np.max(data[:,1], axis=0)
   
    return (mean, low, high)

def accumulate_files(n):
    meanarr = []
    lowarr = [] 
    higharr = []

    for m in range(10,101,5):
        # for g=0
        filename = f"est_tries/data{m}x{n}"
        (mean, low, high) = process_data(filename)
        meanarr += [mean]
        lowarr += [low]
        higharr += [high]

    meanarr = np.array(meanarr)
    lowarr = np.array(lowarr)
    higharr = np.array(higharr)
    mvalues = np.arange(10,101,5) 

    return (meanarr, lowarr, higharr, mvalues)

def plot_tries(mean, low, high, xaxis):

    fig, ax = plt.subplots(1)
    ax.plot(xaxis, mean, lw=1, label="mean", color='black', ls='--')
    ax.fill_between(xaxis, low, high, facecolor='blue', alpha=0.3,
                label='range')

    ax.legend(loc="center right")
    ax.set_xlabel("m (number of rows in P)")
    ax.set_ylabel("estimated # strings to check")
    ax.grid()
    ax.set_yscale('log',basey=2)

    plt.title("Number of Candidate strings to try when n=50")
    plt.show()



# for n = 50
(mean, low, high, xaxis) = accumulate_files(50)

#plot_occ(occ,xaxis)
plot_tries(mean, low, high, xaxis)
