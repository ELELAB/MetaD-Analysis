#!/usr/bin/env python
__author__ = "Mads Nygaard"
__date__ = "20171011"

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import re
import subprocess
import os
from matplotlib.widgets import Slider
import mmap
from contextlib import closing
from itertools import combinations
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.axes_grid1 import make_axes_locatable




def yes_no_input(string): #Is not enabled at the moment. Get yes/no from user.
    is_valid = 0
    while not is_valid:    
         choice = raw_input(string+" [y/n]: ")
         if choice in ["Y","y"]:
             is_valid = 1
             return True
         elif choice in ["N","n"]:
             is_valid = 1
             return False
         else:
             print "Not valid input: " + choice  

#def lastline(fname):
#    with open(fname, 'r+b') as myfile:
#            with closing(mmap.mmap(myfile.fileno(), 0, access=mmap.ACCESS_WRITE)) as mm:
#                startofline = mm.rfind(b'\n', 0, len(mm) - 1) + 1
#                return mm[startofline:].rstrip(b'\r\n')


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')): #A key in order to sort numbers naturally
    return [int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s)]   

def custom_plot(ax, x, y, y2,fill=True,fillC=None,**kwargs): #Plots data to ax, fills below if fill=True
#   ax = kwargs.pop('ax', plt.gca())
    base_line, = ax.plot(x, y, **kwargs)
    if fillC == None:
        fillC = base_line.get_color()
    if fill: ax.fill_between(x,y,y2,where=y>y2[0],facecolor=fillC) #Fills with same color as the line

#A fileparser class with a lot of static methods
class FileParser():
    def __init__(self):
        pass
    def get_data(fesfile):
        pass

    @staticmethod
    def lastline(fname):
        """Returns last line of file as string"""
        with open(fname, 'r+b') as myfile:  #Apparently fast code for getting last line in long files
                with closing(mmap.mmap(myfile.fileno(), 0, access=mmap.ACCESS_WRITE)) as mm:
                    startofline = mm.rfind(b'\n', 0, len(mm) - 1) + 1
                    return mm[startofline:].rstrip(b'\r\n')

    @staticmethod
    def getinfo(fname):
        """Gets info from comment lines in Plumed2 files in dict
        {"FIELDS":[...],"set1":["True/False"]...}"""
        with open(fname,"r") as f:
            trig = True
            comments ={}
            while trig == True: #Reads until no more lines with #!
                line = f.readline().split()
                if "#!" in line:
                    if line[1]=="SET":
                        comments[line[2]]=line[3]
                    else:
                        comments[line[1]]=line[2:]
                else:
                    trig = False
        return comments

    @staticmethod
    def _findminmax(CVfile):
        """Returns min and max X and Y val from fes file
        (min(x),max(x)),(min(y),max(y))"""
        load=np.loadtxt(CVfile)
        x=load[:,0]
        y=load[:,1]
        return (min(x),max(x)),(min(y),max(y))
    pass

class BashRunner(): #Runs the different tools needed for analysis in Bash
    def __init__(self,fn_hills="./HILLS"):
        self.fn_hills = fn_hills
        self.hills = FileParser.getinfo(fname=fn_hills)
        self.CVs = []
        self.simlenght = float(FileParser.lastline(fn_hills).split()[0])
        self.fldn_fes = []
        self.fldn_2d = []
        self.fldn_1d = []
        self.runs = 0
        self.miscplots = []
        for i in self.hills["FIELDS"]:  #Gets the CVs described in the HILLS file
            if "sigma" not in i:
                if "time" not in i:
                    self.CVs.append(i)
            else:
                break

    def multithreader(self,l_processes,n_threads): #A default mulitrheader
        import Queue
        import threading
        import time
        import subprocess
        exitFlag = 0
        class updatedThread(threading.Thread):
            def __init__(self, threadID, name, q):
                threading.Thread.__init__(self)
                self.threadID = threadID
                self.name = name
                self.q = q
            def run(self):
                print "Starting T" + self.name
                process_data(self.name, self.q)
                print "Exiting T" + self.name

        def process_data(threadName, q):
            while not exitFlag:
                queueLock.acquire()
                if not workQueue.empty():
                    data = q.get()
                    queueLock.release()
                    print "T%s running:%s" % (threadName, data)  
                    if data == "S":
                        time.sleep(4)
                    else:
                        time.sleep(1)
                    process = subprocess.Popen(data,stdout=subprocess.PIPE)
                    process.communicate()
                else:
                    pass
                    queueLock.release()

        threadList = range(n_threads)
        nameList = l_processes
        queueLock = threading.Lock()
        workQueue = Queue.Queue(0)
        threads = []
        threadID = 1
        tnow = time.time()
        
        # Create new threads
        for tName in threadList:
            thread = updatedThread(threadID, tName, workQueue)
            thread.start()
            threads.append(thread)
            threadID += 1

        # Fill the queue
        queueLock.acquire()
        for word in nameList:
            workQueue.put(word)
        queueLock.release()

        # Wait for queue to empty
        while not workQueue.empty():
            pass

        # Notify threads it's time to exit
        exitFlag = 1

        # Wait for all threads to complete
        for t in threads:
            t.join()
        print "Exiting Main Thread after %d sec" % (time.time()-tnow)

    def _bash_prep(self,cvs=None,plumedloc="/usr/local/plumed-2.2.0/bin/plumed"
            ,fn_hills="./HILLS",stride=None,out_folder="./2Dfes",tempK = 298.0):
        """Prepares the bash command in order to run sum_hills, and creates the folders for files"""
        if cvs is not None:
            cvs = ",".join(cvs)
        sstride = None if stride is None else "--stride"
        fname = "/fes.dat" if stride is None else "/fes_"
        idw = None if cvs is None else "--idw"
        kT = tempK*0.0083144598
        runlist = [plumedloc,"sum_hills","--hills", self.fn_hills, "--kt", "%.2f"%(kT), idw, cvs, sstride, stride, "--mintozero", "--outfile", out_folder+fname]
        runlist2 = [str(x) for x in runlist if x is not None]
#        print "Running: " + " ".join(runlist2)
        #try:
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)
#        process = subprocess.Popen(runlist2,stdout=subprocess.PIPE)
#        stdout,stderr = process.communicate()
        return runlist2
        
    
    def sum_hills(self,CVs=None,plumedloc="/usr/local/plumed-2.2.0/bin/plumed",stride=None,dryrun = False,nt=1,tempK = 298,no2D=False):
        """Prepares all runs specified by the CVs and runs them with the multithreader."""
        if CVs is None:
            CVs = self.CVs
        if no2D == True:
            self.doublefes = []
        else:
            self.doublefes = list(combinations(CVs,2))
        listofruns=[]
        
 #       if (len(self.doublefes)+len(CVs)) > 10:
 #           if yes_no_input(("Attention, this will run %s different calculations"
 #                   "are you shure you want to continue?")%(len(self.doublefes)+len(CVs))) == False:
 #                   exit()
        if len(self.doublefes) > 1:
            for pair in self.doublefes:
                folder = "./2D_"+pair[0] + "_" + pair[1]
                self.fldn_fes.append(folder)
                self.fldn_2d.append(folder)
                self.runs +=1
                listofruns.append(self._bash_prep(cvs=pair, plumedloc=plumedloc, fn_hills=self.fn_hills, stride=stride, out_folder=folder, tempK = tempK))
#               process = subprocess.Popen([plumedloc,"sum_hills","--hills", self.fn_hills, "--kt", "2.5", "--stride", "20000", "--outfile", folder+"/2Dfes_"],stdout=subprocess.PIPE)
        elif len(self.doublefes) == 1:
            folder = "./2D_"+self.doublefes[0][0] + "_" + self.doublefes[0][1]
            self.fldn_fes.append(folder)
            self.fldn_2d.append(folder)
            self.runs +=1
            listofruns.append(self._bash_prep(cvs=None,plumedloc=plumedloc,fn_hills=self.fn_hills,stride=stride,out_folder=folder,tempK=tempK))
        for CV in CVs:
            folder = "./1D_"+CV
            self.fldn_fes.append(folder)
            self.fldn_1d.append(folder)
            self.runs +=1
            if len(CVs) == 1:
                CV = None
            else:
                CV = [CV]
            listofruns.append(self._bash_prep(cvs=CV,plumedloc=plumedloc,fn_hills=self.fn_hills,stride=stride,out_folder=folder,tempK = tempK))
        #print listofruns
        if dryrun == False:
            self.multithreader(listofruns,nt) 

    def gromacs(self):
        """TODO: Returns Rg and Rh value for plotting, tests for PBC"""
        self.runs+=1
    
    def crysol():
        """TODO: Returns avg Rg"""
        pass

    def fileplotter(self,fn,cols,title=""):
        """Colvar or HILLS file plotter. Specify cols with list of column headings"""
        if not isinstance(cols, list):
            cols = [cols]
        fieldinfo = FileParser.getinfo(fn)["FIELDS"] #Reads the info of columns in file
        load = np.loadtxt(fn)
        t = load[:,0]/1000 
        for col in cols:
            try:
                column = fieldinfo.index(col) #Gets column number
                y = load[:,column]
                self.miscplots.append((np.array((t,y)),("Time[ns]",col),title+" "+col))
                self.runs +=1    
            except KeyError as e:
                print "Column %s is not in %s" %(e.args[0],fn)
        
    def rad_gyr_CA(self,gro,xtc,sel="name CA",stride=1):
        """Plots CA Rg using MDAnalysis, slower than just plotting Colvar.
        Needs the gro and xtc file"""
        try:
            import MDAnalysis
            try:
                u = MDAnalysis.Universe(gro,xtc) #Reads the traj to a universe object 
                time = []
                radius = []
                bb = u.select_atoms(sel) #Select specific atoms
                for ts in u.trajectory[::stride]: #Loads every "stride" frame of the traj
                    time.append(ts.time/1000) #Gets the time of the frame in ns
                    radius.append(bb.radius_of_gyration()) #Calculates the Rg
                
                #np.savetxt("radius.txt", np.array((time,radius)),header="#! time rg_ca")
                toplot = (np.array((time,radius)),("Time[ns]",u"Rg[\xc3]"),"Rg("+sel+")") #Plot is added as touple (data,label,title)
                self.runs +=1 #Adds to runs
                self.miscplots.append(toplot) #Add to plotlist 
                return toplot
                
            except ValueError as e:
                print "Error in Rg calculation:"
                print e
            except TypeError as e:
                print "Something wrong with the xtc or gro filenames:"
                print e
        except ImportError as e:
            print "You need MDAnalysis to run Rgplot, to install run 'pip install --upgrade MDAnalysis'."
            print e
    

class Plotter(): #Used to plot...
    multiPagePlot = False  
    def __init__(self,ntotplot,**kwargs):
        self.ntotplot = ntotplot
        self.sliders = []
        self.contfsliders = []
        self.imax = []
        self.contmeta = []
        self.imlists = []
        self.contlist = []
        self.contax = []
        self.axtime = []
        self.n2d = 0
        self.ncont = 0
        self.fig = plt.figure(facecolor="white")
        self.plotno= 1
        if "nrows" and "ncols" in kwargs:  #Setting up the plot grid
            self.nrows = kwargs["nrows"]
            self.ncols = kwargs["ncols"]
            #self.ax =self.fig.add_subplot(self.nrows,self.ncols,1)    
        else:
            if ntotplot == 1:
                self.ncols = 1
            elif ntotplot <= 6:
                self.ncols = 2
            else:
                self.ncols = 3
            self.nrows = int(ntotplot/self.ncols)+ (1 if ntotplot%self.ncols>0 else 0)

    def _axadder(self,title = None, xylabel = None,**kwargs):
        """Adds an ax to the plot"""
        ax = self.fig.add_subplot(self.nrows,self.ncols,self.plotno,**kwargs)
        ax.tick_params(direction = "out",right=False,top=False)
        self.plotno+=1
        if title is not None:
            ax.set_title(title)
        if xylabel is not None:
            ax.set_xlabel(xylabel[0])
            ax.set_ylabel(xylabel[1])
        return ax

    def simpleplot(self,x,y,title=None,xylabel=None,**kwargs):
        """Just for adding a very simple plot and adding it to plots available"""
        ax = self._axadder(title,xylabel)
        custom_plot(ax,x,y,None,fill=False,**kwargs)
    
    def twoD3Dplot(self, fn_fes, remove=0,rem_dir='x',xlabel=None,ylabel=None,**kwargs):
        """A 3D plotter for 2D fes plotting. Messes up the 'tight_layout' option at the moment.
        Use as a a single plot is best. Cut away 'remove' number of data points from the 'rem_dir'
        direction"""
        from mpl_toolkits.mplot3d import axes3d
        comments = FileParser.getinfo(fn_fes)
        CVs = comments["FIELDS"][0:2]
        ax= self._axadder(title = CVs[0]+" "+CVs[1],xylabel = CVs, projection="3d")
        loadfile = np.loadtxt(fn_fes).T
        xbin = int(comments["nbins_"+CVs[0]]) #381
        ybin = int(comments["nbins_"+CVs[1]]) #53
        x=loadfile[0][:xbin]
        xr = np.tile(x,(ybin,1))
        y=loadfile[1][::xbin]
        yr = np.tile(y,(xbin,1)).T
        Z= loadfile[2].reshape(ybin,xbin) 
        if remove > 0:
            if rem_dir.lower() == 'x':
                newxtrans= xr.T[:remove]
                newytrans= yr.T[:remove]
                newztrans= Z.T[:remove]
                surfx = newxtrans.T
                surfy = newytrans.T
                surfz = newztrans.T
            if rem_dir.lower() == 'y':
                surfx= xr[remove:]
                surfy= yr[remove:]
                surfz= Z[remove:]
        else:
            surfx = xr
            surfy = yr
            surfz = Z
        ax.plot_surface(surfx, surfy, surfz, rstride=3, cstride=16, alpha=0.8,cmap=mpl.cm.jet_r,**kwargs)
        cset = ax.contour(xr, yr, Z,11, zdir='z', offset=np.amin(Z), cmap=mpl.cm.jet_r)
        # Get rid of the panes
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

        # Get rid of the spines
        ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

        # Get rid of the ticks
        ax.grid(False)
        ax.set_zticks([])

        # Add the labels
        ax.set_xlabel(xlabel) if xlabel is not None else ax.set_xlabel(CVs[0])  
        ax.set_ylabel(ylabel) if ylabel is not None else ax.set_ylabel(CVs[1])  
        ax.set_xlim(min(x),max(x))
        ax.set_ylim(min(y),max(y))
        ax.set_zlim(np.amin(Z), np.amax(Z))

    def twoDcontourplot(self,feslist,CVs=None,zeroset = True,**kwargs):
        """Adds a twoD plot with sliders for time to the fig"""
        sortlist = sorted(feslist,key=natural_sort_key)       
        comments = FileParser.getinfo(sortlist[-1])
        if CVs == None:
            CVs = comments["FIELDS"][0:2]
        self.contlist.append([])
        ax= self._axadder(title = CVs[0]+" "+CVs[1],xylabel = CVs)
        
        def update(val):
            for idx,sl in enumerate(self.contfsliders):
                for coll in self.contax[idx].collections:
                    coll.remove()
                
                xr,yr = self.contmeta[idx][0]
                Zmin,Zmax = self.contmeta[idx][1]
                Zdata =self.contlist[idx][int(round(sl.val))]
                self.contax[idx]=self.contax[idx].ax.contour(xr,yr,Zdata,11,vmin=Zmin,vmax=Zmax,cmap=mpl.cm.jet_r,**kwargs)
        
        divider = make_axes_locatable(plt.gca())
        sloc = divider.append_axes("top","5%",pad="1%")
        self.axtime.append(self.fig.add_axes(sloc, axisbg='lightgoldenrodyellow'))
        self.contfsliders.append(Slider(self.axtime[self.ncont],'t', 0, len(feslist)-1
            , valinit=len(feslist)-1,valfmt='%1.0f'))      
        
        xbin = int(comments["nbins_"+CVs[0]]) #381
        ybin = int(comments["nbins_"+CVs[1]]) #53
        lfile = np.loadtxt(sortlist[-1]).T
        x=lfile[0][:xbin]
        y=lfile[1][::xbin]
        xr = np.tile(x,(ybin,1))
        yr = np.tile(y,(xbin,1)).T
        Zmin,Zmax = np.amin(lfile[2]),np.amax(lfile[2])
        Zsub=0
        if zeroset == True:
            Zsub = Zmin
        rZmin,rZmax = Zmin-Zsub,Zmax-Zsub
        self.contmeta.append(([xr,yr],[rZmin,rZmax]))
       
        levels = np.linspace(rZmin,rZmax,num=11)
        for f in sortlist:
            loadfile = np.loadtxt(f).T
            Z= np.subtract(loadfile[2].reshape(ybin,xbin),Zsub)
            self.contlist[self.ncont].append(Z)
        
        self.contax.append(ax.contour(xr, yr, Z,11,vmin=rZmin,vmax=rZmax, offset=np.amin(Z), cmap=mpl.cm.jet_r,**kwargs))
        #self.contax[idx].ax.clabel(self.contax[idx], inline=1, fontsize=10)
        
        cax = divider.append_axes("right","5%",pad = "3%")
        self.fig.colorbar(self.contax[-1],cax=cax)
        self.ncont += 1 
        for s in self.contfsliders:
            s.on_changed(update)
        
    def twoDfesplot(self,feslist,CVs=None): #Old
        """Adds a twoD plot with sliders for time to the fig"""
        sortlist = sorted(feslist,key=natural_sort_key)
        comments = FileParser.getinfo(sortlist[-1])
        if CVs == None:
            CVs = comments["FIELDS"][0:2]
        self.imlists.append([])
        ax = self._axadder(title = CVs[0]+" "+CVs[1],xylabel = CVs)
        def update(val):
            for idx,sl in enumerate(self.sliders):
                self.imax[idx].set_data(self.imlists[idx][int(round(sl.val))])
            #self.ax.canvas.draw()
        #Make the slider:
        divider = make_axes_locatable(plt.gca())
        sloc = divider.append_axes("top","5%",pad="1%")
        self.axtime.append(self.fig.add_axes(sloc, axisbg='lightgoldenrodyellow'))
        self.sliders.append(Slider(self.axtime[self.n2d],'t', 0, len(feslist)-1
            , valinit=len(feslist)-1,valfmt='%1.0f'))      

        for f in sortlist:
            loadfile = np.loadtxt(f).T
            fes_data = loadfile[2].reshape((int(comments["nbins_"+CVs[1]]),int(comments["nbins_"+CVs[0]])))
            self.imlists[self.n2d].append(fes_data)
        
        self.imax.append(ax.imshow(self.imlists[self.n2d][-1]
                ,extent=([float(i) for i in [comments["min_"+CVs[0]]
                ,comments["max_"+CVs[0]],comments["min_"+CVs[1]]
                ,comments["max_"+CVs[1]]]]),aspect="auto"
                ,origin="lower"))
        
        #Make colorbar 
        cax = divider.append_axes("right","5%",pad = "3%")
        self.fig.colorbar(self.imax[-1],cax=cax)

        self.n2d += 1 
        for s in self.sliders:
            s.on_changed(update)

    def oneDfesplot(self,feslist,CV,time=None,zeroset=True,**kwargs):
        """Adds a plot for a list of 1D fes"""
        ax = self._axadder(title=CV,xylabel = (CV,"energy")) 
        
        sortlist = sorted(feslist,key=natural_sort_key)
        colormap = plt.cm.gist_stern_r
        divider = make_axes_locatable(plt.gca())
        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0.0, 1, len(feslist))])
        xval,yval = FileParser._findminmax(sortlist[-1])
        if zeroset == True:
            Ysub = yval[0]
            fillT = True
        else:
            Ysub = 0.0
            fillT = False
        ax.axis(xval+(yval[0]-Ysub,yval[1]-Ysub))
        
        last = len(sortlist)-1
        fillC = None
        for idx,f in enumerate(sortlist):
            if  last == idx:
                fillC = "White"
            load=np.loadtxt(f)
            x=load[:,0]
            y=np.subtract(load[:,1],Ysub)
            y2=np.full(len(y),0*1.05)
            custom_plot(ax,x,y,y2,fill = fillT,fillC = fillC,**kwargs)
        
        cax = divider.append_axes("right","5%",pad="3%")
        timemax = 1 if time is None else time/1000
        norm = mpl.colors.Normalize(vmin=0,vmax=timemax) 
        cb1 = mpl.colorbar.ColorbarBase(cax,cmap=colormap,norm=norm,orientation="vertical")#,ticks=[20,50,60])
        cb1.locator = mpl.ticker.MaxNLocator(nbins = 5)
        cb1.update_ticks()
        cb1.set_label("Time")

    def writePdf(self,fn_pdf,**kwargs):
        """Removes sliders and writes pdf, adds them afterwards again"""
        for slider in self.axtime:
            self.fig.delaxes(slider)
        self.fig.tight_layout(**kwargs)
        print "Writing .pdf"
        self.fig.savefig(fn_pdf)
        for slider in self.axtime:
            self.fig.add_axes(slider)
    
    def multiPage(self,fn_multipage="multiPage.pdf",**kwargs):      
        """Should be possible to plot multiple figures on multiple pages in same pdf, not really working"""
        self.fig.tight_layout(**kwargs)
        if Plotter.multiPagePlot is False:
            Plotter.multiPagePlot = PdfPages(fn_multipage)
        for slider in self.axtime:
            self.fig.delaxes(slider)
        Plotter.multiPagePlot.savefig()
        for slider in self.axtime:
            self.fig.add_axes(slider)
     
    def show(self):
        """Plots the figure in a window"""
        self.fig.tight_layout()
        self.fig.show()      
        raw_input("Press any key to close window")
    pass
"""
bRun= BashRunner(fn_hills="./HILLS")  #Initiate the Bash runner with the Hills file
some = Plotter(2)  #Initiates plotter with number of ax = bRun.runs
bRun.sum_hills(stride="5000",CVs=None,dryrun = "True" )

some.twoDcontourfplot(glob.glob("./2D_ab1_rg_ca"+"/fes*.dat"),bRun.doublefes[0])
some.twoDcontourfplot(glob.glob("./2D_ab1_rg_ca"+"/fes*.dat"),bRun.doublefes[0],zeroset=False)
some.show()
"""
if __name__ == "__main__":
    #Creates parser for user input:
    import argparse 
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description = ("""This is a script that can assist in plotting data of MetaD simulations run with plumed2 with matplotlib. The script runs sum_hills for specified CVs and the the combination pairs of those."""))
    parser.add_argument("-n", action = "store_true", help = "Dry run, no sum_hills calculation. Use if calculation is already run before")
    parser.add_argument("-cvs",nargs="+", metavar = "rg ab", type=str,help= "Specified CVs. If no CVs specified: Sumhills for all CVs+pairs of those will be run.")
    parser.add_argument("-stride",metavar="t [fs]", type = int, help="Stride length")
    parser.add_argument("-temp",metavar="298 [K]",default=298,type=float,help="Specify temperature of simulation in Kelvin")
    parser.add_argument("-hills", metavar="./HILLS",type = str, default="./HILLS", help="Location of Hills file")
    parser.add_argument("-colvar", metavar="./colvar", type = str, default="./colvar", help="Location of the colvar file")
    parser.add_argument("-plumed",metavar="/usr/local/plumed-2.2.0/bin/plumed",type = str, help ="plumed location",default="/usr/local/plumed-2.2.0/bin/plumed")
    parser.add_argument("-plot",action = "store_true", help = "Plot in window")
    parser.add_argument("-pdf",metavar="filename.pdf", type = str, help ="Write .pdf at specified file location")
    parser.add_argument("-nocolvar",action="store_true",help = "Do not plot the CVs from colvar")
    parser.add_argument("-nohills",action="store_true", help = "Do not plot the hillsheight")
    parser.add_argument("-nt",metavar="N", default=1,type = int, help="specify number of threads")
    parser.add_argument("-no2D",action = "store_true",help = "Do not run 2D sum_hills")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-rg",action = "store_true", help = "Plot Rg from .xtc and .gro file")
    group.add_argument("-xtc",metavar="fname.xtc", type = str, help="specify .xtc file")
    group.add_argument("-gro",metavar="fname.gro", type = str, help="specify .gro file")
    args=parser.parse_args()

    #args.no2D = True
    #args.stride=2000
    #args.n=True #Remove this line 
    #args.plot = True #And this one
    #args.cvs = ["rg_ca"]
    #args.rg=False #And this
    #args.xtc="./trj_pcbfix.xtc"
    #args.gro="./protein.gro"
    #args.nohills = True """   
    bRun= BashRunner(fn_hills=args.hills)  #Initiate the Bash runner with the Hills file

    #Test the user input and environmnent before running sum_hills

    if os.path.isfile(args.plumed) == False and args.n == False: #Does plumed exist here?
        raise OSError("Plumed not found on location: %s"%args.plumed)
   

    if (args.plot == True) or (args.pdf is not None): #Is user running with a display?
        try:
            os.environ['DISPLAY']
        except KeyError:
            print "ERROR: No display found, if using ssh please run with -X flag"
            exit(1)
   

    if args.cvs is not None:   #If user has specified CVs, test if they exist
        for cv in args.cvs:
            if cv not in bRun.CVs:
                raise KeyError("The CV:'%s' is not forund in %s. Only the args %s were found"%(cv,args.hills,", ".join(bRun.CVs)))
    
    if args.nocolvar == False: #CVs in colvar, added to plotlist
        bRun.fileplotter(fn=args.colvar,cols=bRun.CVs,title="Colvar")

    #Run the sum hills, if args are None, sum_hills choses defaults. FES plots are added to the plotlist
    bRun.sum_hills(stride=args.stride,CVs=args.cvs,plumedloc=args.plumed, dryrun = args.n, nt=args.nt,tempK=args.temp,no2D = args.no2D)
    

    if args.rg == True: #Runs MDAnalysis Rg calculation, adds to plotlist, saves data in bashrunner
        bRun.rad_gyr_CA(args.gro,args.xtc,stride=10)

     
    if args.nohills == False: #Hillsheight added to plotlist
        bRun.fileplotter(fn=bRun.fn_hills,cols=["height"],title="Hillsheight")

    if (args.plot is True) or (args.pdf is not None): #Actual plotting with matplotlib
	if(bRun.runs <= 6):
            some = Plotter(bRun.runs)  #Initiates plotter with number of ax = bRun.runs
            for plots in bRun.miscplots: #Plots the rg,colvar,nohills
                some.simpleplot(plots[0][0],plots[0][1],title=plots[2],xylabel=plots[1])
            for fld in bRun.fldn_1d: #Plots all 1D FES
                print "Plotting ", fld
                print 
                if glob.glob(fld+"/fes*.dat") == []:
                    raise Exception("No fes files found in %s. Try running with sum_hills again" % fld)                   
                some.oneDfesplot(glob.glob(fld+"/fes*.dat"),fld[5:], time=bRun.simlenght)
            for idx,fld in enumerate(bRun.fldn_2d): #Plots all 2D FES
                print "Plotting ", fld
                if glob.glob(fld+"/fes*.dat") == []:
                    raise Exception("No fes files found in %s. Try running with sum_hills again" % fld)                   
                some.twoDcontourplot(glob.glob(fld+"/fes*.dat"),bRun.doublefes[idx])
        else:
            some = Plotter(len(bRun.miscplots))  #Initiates plotter with number of ax = bRun.runs
            for plots in bRun.miscplots: #Plots the rg,colvar,nohills
                some.simpleplot(plots[0][0],plots[0][1],title=plots[2],xylabel=plots[1])
            some2 = Plotter(len(bRun.fldn_1d))  #Initiates plotter with number of ax = bRun.runs
            for fld in bRun.fldn_1d: #Plots all 1D FES
                if glob.glob(fld+"/fes*.dat") == []:
                    raise Exception("No fes files found in %s. Try running with sum_hills again" % fld)                   
                some2.oneDfesplot(glob.glob(fld+"/fes*.dat"),fld[5:], time=bRun.simlenght)
            some3 = Plotter(len(bRun.fldn_2d))  #Initiates plotter with number of ax = bRun.runs
            for idx,fld in enumerate(bRun.fldn_2d): #Plots all 2D FES
                print "Plotting ", fld
                if glob.glob(fld+"/fes*.dat") == []:
                    raise Exception("No fes files found in %s. Try running with sum_hills again" % fld)                   
                some3.twoDcontourplot(glob.glob(fld+"/fes*.dat"),bRun.doublefes[idx])


#    some.twoD3Dplot("./2D_ab1_rg_ca/fes_5.dat")

    if args.plot is True: #Show the plot in the display
        some.show()
        some2.show()
        some3.show()
    
    if args.pdf is not None: #Print as a pdf
        if(args.pdf.rfind(".") != -1):
            pdf = args.pdf[:args.pdf.rfind(".")]
            ext = args.pdf[args.pdf.rfind("."):]
        else:
            pdf = args.pdf
            ext = ""
        some.writePdf(fn_pdf=pdf+ext)
        some2.writePdf(fn_pdf=pdf+"_1D"+ext)
        some3.writePdf(fn_pdf=pdf+"_2D"+ext)

