
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
from scipy.stats import binned_statistic as bn
from scipy.optimize import curve_fit as cf
from IPython.display import clear_output
import numpy as np
import os
from matplotlib.backends.backend_pdf import PdfPages
from PyAstronomy.pyasl import foldAt
from astropy.stats import LombScargle
from astropy import units as u
import shutil
get_ipython().run_line_magic('matplotlib', 'inline')
o=np.linspace(0,1,10000)
path = '../astars/'
path2 = path + 'Good/'
def curve(theta,a,phi):
    return a*np.sin(2*np.pi*theta-(phi))+np.mean(cv)


# In[ ]:


for root, dirs, files in os.walk(path): #Walks through all the stars in a project directory
    for file in files:
        if (file == "comp.dat"):
            ct = np.loadtxt(root+'/'+file)[0]#Loads composite light curve information
            cv = np.loadtxt(root+'/'+file)[1]
            t0 = np.min(ct)
            if not os.path.isfile(root+'/pgram.dat'): #Computes Generalised L-S periodogram if it doesn't exist
                frequency, power = LombScargle(ct*u.day, cv*u.mag).autopower(minimum_frequency=2/((np.max(ct)-np.min(ct))*u.day))
                np.savetxt(root + '/pgram.dat',(power,frequency))
            else: #If pgram.dat exists, this star has been looked at so skip
                continue
            P=foldAt(ct*u.day,1/frequency[np.argmax(power)]) #Phase-folds on day at maximum power
            maxp = 1/frequency[np.argmax(power)]
            pout,pcov=cf(curve,P,cv) #Least-squares fit sine wave to phase-folded LC
            perr = np.sqrt(np.diag(pcov)) #Square root of diag(covariance) for errors
            bins = bn(P,cv,statistic = 'median',bins = 15) #Sets up 15 bins to plot general form of pfolded LC
            binsites = np.array([])
            for i in range(len(bins[1]) - 1):
                binsites = np.append(binsites, np.median([bins[1][i],bins[1][i+1]]))
            f, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2, 2, sharex=False,sharey = False,squeeze = True,figsize=(10,10))
            
            ax1.set_ylabel("Brightness (mag)") #plots LC
            ax1.set_xlabel("BJD - " + str(t0))
            ax1.set_title("ASCC " + str(regex.findall(root)[0]))
            ax1.invert_yaxis()
            ax1.plot(ct-t0,cv,'black',marker='.',linestyle='None')

            ax2.set_ylabel("Brightness (mag)") #Plots phase folded lc
            ax2.set_xlabel("Phase")
            ax2.set_title(str(np.around(maxp,4)) + " days, " + str(np.around(1/maxp,4)) + " c/d, " +str(np.around(np.abs(pout[0]),4)) +" mag")
            ax2.invert_yaxis()
            ax2.plot(P,cv,'gray',marker='.',linestyle='None')
            ax2.plot(o,curve(o,*pout),'black',linewidth = 5.0)
            ax2.plot(binsites,bins[0],'black',marker='o',markersize=7.0,linestyle='None')
            ax2.set_xlim(0,1)
            #ax1[0].axes.get_xaxis().set_visible(False)

            ax3.set_ylabel("Normalized Power") #Plots time periodogram
            ax3.set_xlabel("Time (Days)")
            ax3.plot(1/frequency,power,'black')
            ax3.tick_params(axis='both', which='minor')

            ax4.set_ylabel("Normalized Power") #Plots frequency periodogram
            ax4.set_xlabel("Frequency (c/d)")
            ax4.plot(frequency,power,'black')
            ax4.tick_params(axis='both', which='minor')
            plt.show()
            
            name = input("Good (l) or bad (else)?") #Initial sorting. Basically does it have a unique period to
            #the distribution of systematics?
            type(name)
            if(name == 'l'):
                clear_output()
                pp=PdfPages(root + '/plots.pdf') #Saves copy of plots to folder
                pp.savefig(f,bbox_inches='tight')
                pp.close()
                plt.close('all')
                shutil.move(root,"../astars/Good") #Moves folder to step 2 (next code segment)
            else:
                clear_output()#Moves onto the next star
                plt.close("all")


# In[ ]:


for root, dirs, files in os.walk(path): #Individual sorting code
    print(root)
    for file in files:
        name = 'q' #Sets intiial variables for a given star in the "Good" folder of a project. name is what is used
        #to process what to do next.
        d = True
        in6 = 'f'
        while(d == True):
            clear_output()
            plt.close("all")
            if (file == "comp.dat"):
                ct = np.loadtxt(root+'/'+file)[0]
                cv = np.loadtxt(root+'/'+file)[1]
                t0 = np.min(ct)
                power,frequency = np.loadtxt(root+'/pgram.dat')
                if name != 'c':
                    maxp = 1/frequency[np.argmax(power)]
                else:
                    maxp = in5
                P=foldAt(ct*u.day,maxp) #Recomputes the phas-fold as needed
                pout,pcov=cf(curve,P,cv)
                perr = np.sqrt(np.diag(pcov))
                bins = bn(P,cv,statistic = 'median',bins = 40)
                binsites = np.array([])
                for i in range(len(bins[1]) - 1):
                    binsites = np.append(binsites, np.median([bins[1][i],bins[1][i+1]]))
                f, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2, 2, sharex=False,sharey = False,squeeze = True,figsize=(10,10))
                #Draws same plot as above, but with possible new options (below).
                ax1.set_ylabel("Brightness (mag)")
                ax1.set_xlabel("BJD - " + str(t0))
                ax1.set_title("ASCC " + str(regex.findall(root)[0]))
                if name == 'z':
                    ax1.set_xlim(in1,in2)
                ax1.invert_yaxis()
                ax1.plot(ct-t0,cv,'black',marker='.',linestyle='None')

                ax2.set_ylabel("Brightness (mag)")
                ax2.set_xlabel("Phase")
                ax2.set_title(str(np.around(maxp,4)) + " days, " + str(np.around(1/maxp,4)) + " c/d, " +str(np.around(np.abs(pout[0]),4)) +" mag")
                ax2.invert_yaxis()
                ax2.plot(P,cv,'gray',marker='.',linestyle='None')
                if not (in6 == 'y'):
                    ax2.plot(o,curve(o,*pout),'black',linewidth = 5.0)
                ax2.plot(binsites,bins[0],'black',marker='o',markersize=7.0,linestyle='None')
                ax2.set_xlim(0,1)

                ax3.set_ylabel("Normalized Power")
                ax3.set_xlabel("Time (Days)")
                ax3.plot(1/frequency,power,'black')
                if name == 'x':
                    ax3.set_xlim(in3,in4)
                ax3.tick_params(axis='both', which='minor')

                ax4.set_ylabel("Normalized Power")
                ax4.set_xlabel("Frequency (c/d)")
                ax4.plot(frequency,power,'black')
                if name == 'x':
                    if(in3 == 0):
                        ax4.set_xlim(0,1/in4)
                    else:
                        ax4.set_xlim(1/in3,1/in4)
                ax4.tick_params(axis='both', which='minor')
                plt.show()
                if name =='r':
                    print(newmax)
                name = input("z = zoom light curve, x = zoom pgram, c = phase-fold d, s = save plot, v= remove,\n f = cepheid, g = EB, r = search, \np=DS, y = LV, i = rot, u = irregular, j = Be, else=next")
                type(name)
            if(name == 'z'): #Redraw with new LC time limits
                in1 = float(input("Lower limit = "))
                type(in1)
                in2 = float(input("Upper limit = "))
                type(in2)
            elif(name == 'x'): #Redraw pgrams with new time/frequency bounds
                in3 = float(input("Lower limit = "))
                type(in3)
                in4 = float(input("Upper limit = "))
                type(in4)
            elif(name == 'c'): #Phase-fold on new period
                in5 = float(input("New Period = "))
                type(in5)
                in6 = input("Remove Sinusoid (y)?")
                type(in6)
            elif(name == 'r'): #Searches for the stronges period in the given range of days.
                in7 = float(input("Lower limit = "))
                type(in7)
                in8 = input("Upper limt = ")
                type(in8)
                days2 = 1/frequency[np.where(np.logical_and(1/frequency>=float(in7),1/frequency<=float(in8)))]
                power2 = power[np.where(np.logical_and(1/frequency>=float(in7),1/frequency<=float(in8)))]
                newmax = days2[np.argmax(power2)]
            elif(name == 's'): #Saves current plot over plots file in folder
                pp=PdfPages(root + '/plots.pdf')
                pp.savefig(f,bbox_inches='tight')
                pp.close()
            elif(name == 'v'): #Each of these sorts into respective variable type.
                shutil.move(root,path)
                d = False
            elif(name == 'f'):
                shutil.move(root,'../Cepheids')
                d = False
            elif(name == 'y'):
                shutil.move(root,'../Long_Variable')
                d = False
            elif(name == 'p'):
                shutil.move(root,'../DS')
                d = False
            elif(name == 'u'):
                shutil.move(root,'../Irregular')
                d = False
            elif(name == 'i'):
                shutil.move(root,'../Rot')
                d = False
            elif(name == 'g'):
                shutil.move(root,'../EBs')
                d = False
            elif(name == 'j'):
                shutil.move(root,'../BE')
                d = False
            else: 
                d = False

