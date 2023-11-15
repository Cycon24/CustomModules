# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 13:53:19 2023

@author: cycon

"""
import numpy as np
import matplotlib.pyplot as plt

def fft(xlow, xhigh, y_funct, N=300, barwidth=0.2):
    T = 1/N
    
    x = np.linspace(xlow,xhigh,N, endpoint=False)
    
    if ~callable(y_funct):
        print('Error: Inputted Function is not Callable')
        return None
    else:
        y = y_funct(x)
        
    fig, ax = plt.subplots(1,2)
    fig.set_size_inches([9,5])
    ax[0].plot(x,y)
    ax[0].set_xlabel('x')
    ax[0].set_ylabel('y')
    
    amp = abs(np.fft.fft(y))#[:N//2]
    freq = np.fft.fftfreq(N,T)#[:N//2]
    
    ax[1].bar(freq, 2*amp/N, width=barwidth)
    ax[1].plot(freq, 2*amp/N, 'r')
    ax[1].set_xlabel('frequency')
    ax[1].set_ylabel('amplitude')
    ax[1].set_xlim(0,30)
    return

if __name__=="__main__":
    def y(x): 
        return 2*np.sin(2*np.pi*x/7.)-4*np.sin(3.*np.pi*x/5.)
    fft(0,9, y)