import numpy as np
import matplotlib.pyplot as plt


# prg seed, status (0 success, 1 ran out of attempts, -1 P not full rank, 2 wrong string found), #attempts

def process_data(filename):
    data = np.genfromtxt(filename,delimiter=',', dtype=np.float64)
    #strip the prg seed
    data = data[:,1:-1]
    print(data.shape)
    #count how many with each status
    mean = np.mean(data,axis=0)
    low = np.min(data, axis=0)
    high = np.max(data, axis=0)
    
    return (mean, low, high)


def plot_rank(mean, low, high, xaxis):

    fig, ax = plt.subplots(1)
    ax.plot(xaxis, mean, lw=1, label="mean", color='black', ls='--')
    ax.fill_between(xaxis, low, high, facecolor='blue', alpha=0.3,
                label='range')

    ax.legend(loc="center right")
    ax.set_xlabel("k")
    ax.set_ylabel("rank Q")
    ax.grid()
    ax.set_yscale('log',base=2)
    ax.set_xticks(xaxis)
    plt.title("rank Q when n,m = 500")
    plt.show()


filename = f"rankQ/data500x500"
# for n = 50
(mean, low, high) = process_data(filename)
xaxis = np.arange(1,11);

plot_rank(mean, low, high, xaxis)
