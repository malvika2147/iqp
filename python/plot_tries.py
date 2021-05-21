import numpy as np
import matplotlib.pyplot as plt


# prg seed, status (0 success, 1 ran out of attempts, -1 P not full rank, 2 wrong string found), #attempts

def process_data(filename):
    data = np.genfromtxt(filename,delimiter=',', dtype=np.float64)
    #strip the prg seed
    data = data[:,1:]

    #count how many with each status
    succ  = np.count_nonzero(data[:,0] == 0)
    fail =  np.count_nonzero(data[:,0] == 1)
    wrong = np.count_nonzero(data[:,0] == 2)
    skip = np.count_nonzero(data[:,0] == -1)
    #number of tests
    ntests = data.shape[0]; 

    # warn if code ran out of attempts 
    assert fail == 0
    
    # NOTE: fix percentages if ntests changes
    assert ntests == 100 
    # NOTE: order is succ, wrong, skip
    arr = np.array([succ,wrong,skip])
    
    if skip == 100:
        return (arr, -1, -1, -1)

    mean = np.mean(data[:,1], where=data[:,0]>=0,axis=0)
    low = np.min(data[:,1], where=data[:,0]>=0, axis=0, initial = 1e10)
    high = np.max(data[:,1], where=data[:,0]>=0,axis=0, initial = -1)
   
    return (arr, mean, low, high)

def accumulate_files(n):
    occ = []
    meanarr = []
    lowarr = [] 
    higharr = []
    xfil = []

    for m in range(40,301,10):
        # for g=0
        filename = f"tries/data{m}x{n}"
        (counts, mean, low, high) = process_data(filename)
        occ += [counts]

        # ignore tests where everything was skipped
        if mean >= 0:
            meanarr += [mean]
            lowarr += [low]
            higharr += [high]
            xfil += [m]

    meanarr = np.array(meanarr)
    lowarr = np.array(lowarr)
    higharr = np.array(higharr)
    occ = np.array(occ)
    mvalues = np.arange(40,301,10) 
    xfil = np.array(xfil)

    return (occ, meanarr, lowarr, higharr, mvalues, xfil)

def plot_occ(occ, xaxis):
    
    fig, ax = plt.subplots(1)
    ax.plot(xaxis, occ[:,0], 'o',  label="successful", color='green')
    ax.plot(xaxis, occ[:,2], '^', label="skipped", color='black')
    ax.plot(xaxis, occ[:,1], 'x',  label="wrong string output", color='red')
    
    ax.legend(loc="center right")
    ax.set_xlabel("m (number of rows in P)")
    ax.set_ylabel("% of tests")
    ax.grid()

    plt.title("Status of searching for s when n=50")
    plt.show()


def plot_tries(mean, low, high, xaxis):

    fig, ax = plt.subplots(1)
    ax.plot(xaxis, mean, lw=1, label="mean", color='black', ls='--')
    ax.fill_between(xaxis, low, high, facecolor='blue', alpha=0.3,
                label='range')

    ax.legend(loc="center right")
    ax.set_xlabel("m (number of rows in P)")
    ax.set_ylabel("# strings tried per test")
    ax.grid()
    ax.set_yscale('log',basey=2)

    plt.title("Number of Candidate string tried when n=50")
    plt.show()



# for n = 50
(occ, mean, low, high, xaxis, xfilt) = accumulate_files(50)

#plot_occ(occ,xaxis)
print(mean)
plot_tries(mean, low, high, xfilt)
