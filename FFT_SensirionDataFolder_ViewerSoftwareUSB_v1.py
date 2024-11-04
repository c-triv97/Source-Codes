import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import linecache 
import os
from os import listdir
from os.path import isfile, join
 
 ##### IMPORTANT TO NOTE: THERE IS A SUBSETTING METHOD IN THE START WHICH GENERATES A DIFFERENT 
 # RANDOM SUBSET EACH TIME! ANY REMAINING NAN VALUES ARE INTERPOLATED BEFORE FFT TRANSFORMATION - 
 # this is currently commented out 
  
class SensirionExtract:
    def __init__(self, folder_file) -> None:

        '''
        Python thing, set up the file and extract the data needed
        '''
        self.folder_file=folder_file
        
        print(self.folder_file[:-4])

        #read file into Dataframe
        self.df = pd.read_csv(folder_file, skiprows=15, sep=",") # change to "," if recorded on UK laptop
        print(f"file: {self.folder_file}")

        #extract Hz from file
        self.SAMPLERATE = int(float(linecache.getline(self.folder_file, 12).split(",")[-1].split('\n')[0].replace('"', ''))) # change to "," if recorded on UK laptop
        self.SAMPLERATE = 1/(self.SAMPLERATE/1000)
        print(f"Sample Rate: {self.SAMPLERATE}")
        
        # # subset of random 10s
        # self.df['time_change'] = self.df["Relative Time[s]"].diff().round(4) # for some reason python wont recognise the numbers as exactly 0.001 unless you do this 
        # # choose n row that corresponds to 10s sample 
        # secs=30*1000 
        # # setting up filter conditions for start point
        # c1=self.df["Flow linearized [ml/min]"].notna()
        # c2=self.df["time_change"] == 0.0010
        # # select index of a random row meeting both conditions 
        # start_idx = self.df[c1 & c2].sample(n=1).index[0]
        
        # self.df = self.df.iloc[start_idx:start_idx + secs]
        
        # Count NaNs in each column
        self.nan_count = self.df.isna().sum()
        print("NaN count per column:")
        print(self.nan_count)
        
        # setting relative times 
        self.df["Time"] = self.df["Relative Time[s]"]-self.df.iloc[0]["Relative Time[s]"]

        # print(f"File sample rate: {self.SAMPLERATE}") #debugging
        self.datafft=None
        self.datacut=None
        
    def run(self,fCutLow,fCutHigh,factor, plotlimsTime=None, plotlimsFreq=None, plotlimsFreqmag=None,saveplot=None, outputcsv=True):
        '''
        Runs the file though filter and plotter
        '''
        # Filter it
        self.filter(fCutLow=fCutLow, fCutHigh=fCutHigh, factor=factor) 
        #Save data
        if outputcsv:
            
            outputpath="output"
            outfilename=self.folder_file.split("\\")[1][:-4] #the -4 removes .csv from the string
            
            if not os.path.exists(outputpath):
                os.makedirs(outputpath)
            
            dataout = pd.DataFrame(self.df["Time"])
            dataout["Time"] = self.df["Time"]
            dataout["Raw Data"] = self.df["Flow linearized [ml/min]"]
            dataout["Filtered Data"] = self.datacut
            dataout.to_csv(f"{outputpath}/{outfilename}_filtered.csv",index_label = 'Sample Number')

        # Plot it
        self.plot(plotlimsTime=plotlimsTime, plotlimsFreq=plotlimsFreq, plotlimsFreqmag=plotlimsFreqmag,saveplot=saveplot)

    def filter(self,fCutLow, fCutHigh = None, factor =0): 
        '''
        # Sets up the filter and checks paramters are within lims and filterable
        # either pass as single arguments, or as lists i.e for bandstops 0->2 & 5->7 : fCutLow = [0,5], fCutHigh = [2,7], factor =[0,0]
        '''
        # Do the FFT
        self.datafft = np.fft.fft(self.df["Flow linearized [ml/min]"].interpolate()) 
        
        # print(f"fcuthigh: {fCutHigh}") # debugging

        # If a single argument for each is passed just filter once.
        if not isinstance(fCutLow, list) and not isinstance(fCutHigh, list) and not isinstance(factor, list):
            # Check cutoffs are within lims of the sample rate
            fCutHigh = self.SAMPLERATE//2 if ((fCutHigh=="Nq") or (fCutHigh == None)) else fCutHigh
            fCutHigh = self.SAMPLERATE//2 if (fCutHigh > self.SAMPLERATE//2) else fCutHigh
            fCutLow  = self.SAMPLERATE//2 if (fCutLow>self.SAMPLERATE//2) else fCutLow
        
            #Filter the data
            datfftcut = self.freqalter(self.datafft,fCutLow,fCutHigh,factor)

        # If a multiple arguments for each are passed recursivly filter it. - not actually tested yet.
        else:
            #check that correct number of parameters passed
            if isinstance(fCutLow, list) and isinstance(fCutHigh, list) and isinstance(factor, list):
                if len(fCutLow) != len(fCutHigh) != len(factor):
                    raise Exception("Length of passed cut off values do not match")
            else:
                raise Exception("Not all arguments are lists!")
            
            datfftcut=self.datafft
            for low,high,fac in zip(fCutLow,fCutHigh,factor):
                # Check each cutoff is within lims of the sample rate
                high = self.SAMPLERATE//2 if ((high=="Nq") or (high == None)) else high
                high = self.SAMPLERATE//2 if (high > self.SAMPLERATE//2) else high
                low  = self.SAMPLERATE//2 if (low>self.SAMPLERATE//2) else low

                #Filter the data
                datfftcut = self.freqalter(datfftcut,low,high,fac)
                

        #Turn frequency domain signal back to time domain
        self.datacut = np.real(np.fft.ifft(datfftcut))

        return self.datacut

    def freqalter(self, data, fstart, fstop, factor):
        '''
        Alters the correct frequncy bins in the FFT array between fstart -> fstop by factor. 
        Singal is real valued signal, so both base and mirror need changed, to result in real valued signal once inversed
        @param factor - the value to multiple those frequesncies by, i.e 0 = remove, 0.5 = halve signal, 2= double signal
        
        '''
        # print(f"params: {fstart} {fstop} {factor} {self.SAMPLERATE}") # debugging
        
        # Find lower and upper bounds of FFT array that contains the frequencies of interest
        k1 = int(fstart*len(data)/self.SAMPLERATE)
        k2 = int(fstop*len(data)/self.SAMPLERATE)

        # Calculate bandwidth of samples
        freqbw = k2-k1
        datamanip = data[:] # stop it changing the array inplace - silly python

        # Scroll through data array and alter data of relevent frequencies
        for i in range(freqbw):
            # For real freuqnecies
            datamanip[k1+i] = factor*datamanip[k1+i] 
            # For Mirror freuqnecies
            datamanip[len(datamanip)-k2+i] = factor*datamanip[len(datamanip)-k2+i]
        return datamanip
    
    def plot(self, plotlimsTime=None,plotlimsFreq=None,plotlimsFreqmag=None, saveplot=False):
        '''
        Plots the filtered and unfiltered time domain signal and the frequency domain signal of the filtered data
        @param plotlimsTime - [start,stop] - zoom in on section of time data, makes easier to see rises and falls
        @param plotlimsFreq = [startf,stopf] - Frequnecies to zoom in on (best set from 0 to Niquist)
        @param plotlimsFreqmag - value to make FFT plot look nicer (Huge DC spike will make frequncies hard to see)
        @param saveplot - Saves the plot, or not
        '''
        if plotlimsFreq[1] == "Nq" or None:
            plotlimsFreq = [0,self.SAMPLERATE//2]
                
        plt.figure(figsize=(12, 8))
        plt.subplot(3, 1, 1)
        plt.plot(self.df["Time"],self.df["Flow linearized [ml/min]"],color='tab:blue')
        plt.title(f'{self.folder_file}')
        plt.xlabel('Time [s]')
        plt.ylabel('Flow Rate')
        plt.xlim(plotlimsTime)

        plt.subplot(3, 1, 2)
        plt.plot(self.df["Time"], self.datacut, label='Filtered', color='tab:red')
        plt.title('Filtered Data')
        plt.xlabel('Time [s]')
        plt.ylabel('Flow Rate')
        plt.xlim(plotlimsTime)

        datcontfft = np.fft.fft(self.datacut)
        freq_steps1 = np.fft.fftfreq(self.datacut.size, d=1/self.SAMPLERATE)

        plt.subplot(3, 1, 3)
        plt.plot(freq_steps1,np.abs(datcontfft),color='tab:green')
        plt.title('Power Spectrum')
        plt.xlabel('Freq [Hz]')
        plt.ylabel('Magnitude [arb]')
        plt.xlim(plotlimsFreq)
        plt.ylim(plotlimsFreqmag)
        plt.tight_layout()
        if saveplot:
            plt.savefig(f"{self.folder_file[-4]}.png", dpi=1000)
        plt.show()


if __name__ =="__main__":
    import glob

    fCutLow = 2 # Highest freuqnecy of interest - anything between this and FcutHigh will be removed; can also be an array of high i.e. [2, 5]
    fCutHigh = "Nq" # set to any value, or None / "Nq" for setting to nyquist; can also be an array of high i.e. [3, "Nq"] - the example here removes 2-3 in addition to 5-Nq
    factor = 0 # 0 removes the data, <0 dampens and >1 enhances 
    plotlimsTime = [0,20] # lengh of time to plot time domain signals
    plotlimsFreq = [0,"Nq"]  # Fdomain plot lims
    plotlimsFreqmag = [0,0.1e6] # get rid of huge DC spike in the plot
    saveplots = False
    outputcsv = True

    folder = "" # folder containing of all the data

    # list all the files in the filder
    files = glob.glob(f"{folder}/*")

    # Go through folder and run each file.
    for file in files:
        do = SensirionExtract(file)
        do.run(fCutLow=fCutLow, fCutHigh=fCutHigh, factor=factor,plotlimsTime=plotlimsTime,plotlimsFreq=plotlimsFreq,plotlimsFreqmag=plotlimsFreqmag, saveplot=saveplots, outputcsv=outputcsv)
        