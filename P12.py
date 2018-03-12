import numpy as np
import multiprocessing as mp
import time
import matplotlib.pyplot as plt

# Sleep for t seconds
def burnTime(t):
    time.sleep(t)

# Main
if __name__ == '__main__':
    N = 16 # The number of jobs
    P = 4  # The number of processes

    # A thread pool of P processes
    pool = mp.Pool(P)
    
    
    # Use a variety of wait times
    wait_time=([10**(-6),10**(-5),10**(-4),10**(-3),10**(-2),10**(-1),10,100])
    ratio = np.zeros(len(wait_time))
    time_ser=np.zeros(len(wait_time))
    time_para=np.zeros(len(wait_time))
    val=np.zeros(len(wait_time))
    for t in range(0,len(wait_time)):
        val=wait_time[t]*np.ones(16)
        # Compute jobs serially and in parallel
        # Use time.time() to compute the elapsed time for each
        serialTime = 1;
        parallelTime = 1;
        # serial
        t1=time.time()
        for j in range(16):
            burnTime(wait_time[t])
        time_ser[t]=time.time()-t1 
       
        t2=time.time()
        pool.map(burnTime,val)
        time_para[t]=time.time()-t2  
        ratio[t]=time_ser[t]/time_para[t]
        print(wait_time[t])
        # Compute the ratio of these times
        #ratio.append(serialTime/parallelTime)

    # Plot the results
    plt.plot(wait_time, ratio, '-ob')
    plt.xscale('log')
    plt.xlabel('Wait Time (sec)')
    plt.ylabel('Serial Time (sec) / Parallel Time (sec)')
    plt.title('Speedup versus function time')
    plt.show()
