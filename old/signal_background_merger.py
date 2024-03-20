import pyhepmc
from pyhepmc.io import WriterAscii
import argparse

import numpy as np
import pandas as pd

from datetime import datetime
from operator import itemgetter
import time
import sys

import psutil

# =============================================================
class signal_background_merger:
    """
    Combine signal and up to three background HEPMC files.
    
    Typical usage:
    python3 signal_background_merger --signalFile dis.hepmc --signalFreq 0 \
            --bg1File hgas.hepmc --bg1Freq 1852 \
            --bg2File egas.hepmc --bg2Freq 1852 \
            --bg3File sr.hepmc
    """

    # ============================================================================================
    def __init__(self):
        """ Parse arguments, print banner, open files, initialize rng"""
        
        self.digestArgs()
        self.rng=np.random.default_rng(seed=self.args.rngSeed)
        self.banner()

    # ============================================================================================
    def digestArgs(self):
        """Handle the command line tedium"""
        parser = argparse.ArgumentParser(description='Merge signal events with up to three background sources.')
        
        parser.add_argument('-i','--signalFile', default='small_ep_noradcor.10x100_q2_10_100_run001.hepmc',
                            help='Name of the HEPMC file with the signal events')
        parser.add_argument('-sf','--signalFreq', type=float, default=0.0,
                            help='Poisson-mu of the signal frequency in ns. Default is 0 to have exactly one signal event per slice. Set to the estimated DIS rate to randomize.')

    
        parser.add_argument('-bg1','--bg1File', default='small_hgas_100GeV_HiAc_25mrad.Asciiv3.hepmc',
                            help='Name of the first HEPMC file with background events')
        parser.add_argument('-bf1','--bg1Freq', type=float, default=31347.0,
                            help='Poisson-mu of the first background frequency in ns. Default is the estimated hadron gas rate at 10x100. Set to 0 to use the weights in the corresponding input file.')

        parser.add_argument('-bg2','--bg2File', default='small_beam_gas_ep_10GeV_foam_emin10keV_vtx.hepmc',
                            help='Name of the second HEPMC file with background events')
        parser.add_argument('-bf2','--bg2Freq', type=float, default=333.0,
                            help='Poisson-mu of the second background frequency in ns. Default is the estimated electron  gas rate at 10x100. Set to 0 to use the weights in the corresponding input file.')

        parser.add_argument('-bg3','--bg3File', default='small_SR_single_2.5A_10GeV.hepmc',
                            help='Name of the third HEPMC file with background events')
        parser.add_argument('-bf3','--bg3Freq', type=float, default=0,
                            help='Poisson-mu of the third background frequency in ns. Default is 0 to use the weights in the corresponding input file. Set to a value >0 to specify a poisson mu instead.')


        parser.add_argument('-o','--outputFile', default='',
                            help='Specify the output file name. By default it will be auto-generated.')

        parser.add_argument('-w','--intWindow', type=float, default=2000.0,
                            help='Length of the integration window in nanoseconds. Default is 2000.')
        parser.add_argument('-N','--nSlices', type=int, default=-1,
                            help='Number of sampled time slices ("events"). Default is 10. If set to -1, all events in the signal file will be used and background files cycled as needed.')
        parser.add_argument('--squashTime', action='store_true',
                            help='Integration is performed but no time information is associated to vertices.')

        parser.add_argument('--rngSeed', action='store',type=int, default=None,
                            help='Random seed, default is None')

        parser.add_argument('-v','--verbose', action='store_true',
                            help='Display details for every slice.')
        self.args = parser.parse_args()

    # ============================================================================================
    def banner(self):
        """Print a banner."""
        print('==================================================================')
        print('=== EPIC HEPMC MERGER ===')
        print('authors: Benjamen Sterwerf* (bsterwerf@berkeley.edu), Kolja Kauder** (kkauder@bnl.gov), Reynier Cruz-Torres***')
        print('* University of California, Berkeley')
        print('** Brookhaven National Laboratory')
        print('*** formerly Lawrence Berkeley National Laboratory')

        print('\nFor more information, run \npython signal_background_merger.py --help')
        print('==================================================================\n')

        freqTerm = str(self.args.signalFreq) + " ns" if self.args.signalFreq> 0 else "(one event per time slice)"
        print('Signal events file and frequency:')
        print('\t- {} \t {}'.format( self.args.signalFile, freqTerm ))
        
        print('\nBackground files and their respective frequencies:')
        if self.args.bg1File != "" :
            freqTerm = str(self.args.bg1Freq) + " ns" if self.args.bg1Freq> 0 else "(from weights)"
            print('\t- {} \t {}'.format( self.args.bg1File, freqTerm ))
        if self.args.bg2File != "" :
            freqTerm = str(self.args.bg2Freq) + " ns" if self.args.bg2Freq> 0 else "(from weights)"
            print('\t- {} \t {}'.format( self.args.bg2File, freqTerm ))
        if self.args.bg3File != "" :
            freqTerm = str(self.args.bg3Freq) + " ns" if self.args.bg3Freq> 0 else "(from weights)"
            print('\t- {} \t {}'.format( self.args.bg3File, freqTerm ))


        print()
        print('Integration Window:',self.args.intWindow, ' ns')
        print('Number of time slices:', self.args.nSlices if self.args.nSlices>0 else ' (as many as needed to exhaust the signal file)')
        if self.args.squashTime:
            print('No timing information will be attached to vertices')
        print('The RNG is seeded to', self.args.rngSeed)
        print('==================================================================\n')

    # ============================================================================================
    def merge( self ):
        """Main method to steer the merging."""

        t0 = time.time()
        # Open signal file
        try :
            # sigFile=pyhepmc.io.ReaderAscii(self.args.signalFile)
            sigFile=pyhepmc.io.open(self.args.signalFile)
        except IOError as e:
            print ('Opening files failed: %s' % e.strerror)
            sys.exit()

        # Open input file readers, group them with the associated frequency,
        # then sort them into signal, frequency bg, or weight type bg
        self.sigDict=dict()
        self.freqDict=dict()
        self.weightDict=dict()
        self.makeDicts ( self.args.signalFile, self.args.signalFreq, signal=True )
        self.makeDicts ( self.args.bg1File, self.args.bg1Freq )
        self.makeDicts ( self.args.bg2File, self.args.bg2Freq )
        self.makeDicts ( self.args.bg3File, self.args.bg3Freq )

        # Open the output file
        if self.args.outputFile != "" :
            outputFileName = self.args.outputFile
        else :
            outputFileName = self.nameGen()

        t1 = time.time()        
        print('Initiation time:',np.round((t1-t0)/60.,2),'min')
                    
        # Process and write
        nSlices=self.args.nSlices

        print()
        print('==================================================================')
        print()
        i=0
        with WriterAscii(outputFileName) as f:
            while True :
                if (i==nSlices) :
                    print ( "Finished all requested slices." )
                    break
                # if  i % 10000 == 0 : print('Working on slice {}'.format( i+1 ))
                if  i % 100 == 0 or self.args.verbose : self.squawk(i)
                hepSlice = self.mergeSlice( i )
                ### Arrgh, GenEvent==None throws an exception
                try : 
                    if hepSlice==None : 
                        print ( "Exhausted signal source." )
                        break
                except TypeError as e:
                    pass # just need to suppress the error

                hepSlice.event_number=i
                f.write_event(hepSlice)
                i=i+1

        t2 = time.time()
        self.slicesDone=i
        print('Slice loop time:',np.round((t2-t1)/60.,2),'min')
        print(' -- ',np.round((t2-t1) / self.slicesDone ,3),'sec / slice')
        

        # Clean up, close all input files
        File, Freq = self.sigDict[self.args.signalFile]
        # File.close()
        for fileName in self.freqDict :
            File, Freq = self.freqDict[fileName]
            File.close()

    
    # ============================================================================================
    def makeDicts( self, fileName, freq, signal=False ):
        """Create background timeline chunks, open input file, sort into frequency or weight type"""
        if fileName == "" : return

        try :
            # Should use File=pyhepmc.open(fileName) but that doesn't currently work
            File=pyhepmc.open(fileName)
        except IOError as e:
            print ('Opening {} failed: {}', fileName, e.strerror)
            sys.exit()
            
        # For the signal only, we keep frequency even if it's 0
        if signal :
            self.sigDict[fileName] = [ File, freq ]
            return
            
        if freq<=0 :
            # file has its own weights
            # In this case, we will need to read in all events
            # because we need the distribution to draw from them
            print ( "Reading in all events from", fileName )
            # remove events with 0 weight - note that this does change avgRate = <weight> (by a little)
            container = [[event, event.weight()] for event in File if event.weight() > 0]
            # sorting may help for lookup, sampling
            container = sorted ( container, key=itemgetter(1) )            
            File.close()
            events  = np.array([ item[0] for item in container ])
            weights = np.array([ item[1] for item in container ])
            avgRate = np.average (weights) # I _think_ this is the total flux
            avgRate *= 1e-9 # convert to 1/ns == GHz
            print( "Average rate is", avgRate, "GHz")
            # np.Generator.choice expects normalized probabilities
            probs = weights / weights.sum()
            self.weightDict[fileName]=[ events, probs, avgRate ]
            return

        self.freqDict[fileName] = [ File, freq ]

        return

    # ============================================================================================
    def mergeSlice(self, i):
        """Arrange the composition of an individual time slice"""
        
        hepSlice = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM)
        
        # Signal and frequency background are handled very similarly,
        # handle them in one method and just use a flag for what little is special
        try :
            hepSlice = self.addFreqEvents( self.args.signalFile, hepSlice, True )
        except EOFError :
            return
        
        # Treat frequency backgrounds very similarly 
        for fileName in self.freqDict :
            hepSlice = self.addFreqEvents( fileName, hepSlice )

        for fileName in self.weightDict :
            hepSlice = self.addWeightedEvents( fileName, hepSlice )

        return hepSlice

    # ============================================================================================
    def addFreqEvents( self, fileName, hepSlice, signal=False ):
        """Handles the signal as well as frequency-style backgrounds"""

        if signal :
            File, Freq = self.sigDict[fileName]
        else :
            File, Freq = self.freqDict[fileName]


        c_light  = 299.792458 # speed of light = 299.792458 mm/ns to get mm
        squash = self.args.squashTime

        # First, create a timeline
        intTime=self.args.intWindow
        # this could be offset by iSlice * intTime for continuous time throughout the file

        # Signals can be different
        if Freq==0:
            if not signal :
                print( "frequency can't be 0 for background files" )
                sys.exit()
            # exactly one signal event, at an arbigtrary point
            slice = np.array([ self.rng.uniform(low=0, high=intTime) ])
        else:
            # Generate poisson-distributed times to place events
            slice = self.poissonTimes( Freq, intTime )

        if self.args.verbose : print ( "Placing",slice.size,"events from", fileName )
        if slice.size == 0 : return hepSlice

        # Insert events at all specified locations
        for Time in slice :
            ### Arrgh, GenEvent==None throws an exception
            inevt = File.read()
            try : 
                if inevt==None :
                    if signal :
                        # Exhausted signal events
                        raise EOFError
                    else :
                        # background file reached its end, reset to the start
                        print("Cycling back to the start of ", fileName )
                        File.close()
                        # File=pyhepmc.io.ReaderAscii(fileName)
                        File=pyhepmc.io.open(fileName)
                        # also update the dictionary
                        self.freqDict[fileName] = [File, Freq]
                        inevt = File.read()
            except TypeError as e:
                pass # just need to suppress the error

            # Unit conversion
            TimeHepmc = c_light*Time

            particles, vertices = [], []
            # Stores the vertices of the event inside a vertex container. These vertices are in increasing order so we can index them with [abs(vertex_id)-1]
            for vertex in inevt.vertices:
                position=vertex.position
                if not squash:
                    position=position+pyhepmc.FourVector(x=0,y=0,z=0,t=TimeHepmc)
                v1=pyhepmc.GenVertex(position)
                vertices.append(v1)
            
            # copies the particles and attaches them to their corresponding vertices
            for particle in inevt.particles:
                # no copy/clone operator...
                momentum, status, pid = particle.momentum, particle.status, particle.pid                
                p1 = pyhepmc.GenParticle(momentum=momentum, pid=pid, status=status)
                p1.generated_mass = particle.generated_mass
                particles.append(p1)
                
                # since the beam particles do not have a production vertex they cannot be attached to a production vertex
                if particle.production_vertex.id < 0:
                    production_vertex=particle.production_vertex.id
                    vertices[abs(production_vertex)-1].add_particle_out(p1)
                    hepSlice.add_particle(p1)
                    
                # Adds particles with an end vertex to their end vertices
                if particle.end_vertex:
                    end_vertex = particle.end_vertex.id
                    vertices[abs(end_vertex)-1].add_particle_in(p1)

            # Adds the vertices with the attached particles to the event
            for vertex in vertices:
                hepSlice.add_vertex(vertex)

        return hepSlice
        
    # ============================================================================================
    def addWeightedEvents( self, fileName, hepSlice, signal=False ):
        """Handles weighted backgrounds"""

        events, probs, avgRate = self.weightDict[fileName]

        c_light  = 299.792458 # speed of light = 299.792458 mm/ns to get mm
        squash  = self.args.squashTime
        intTime = self.args.intWindow
        squash = self.args.squashTime

        # How many events? Assume Poisson distribution
        while True:
            nEvents = int( self.rng.exponential( intTime * avgRate ) )
            if nEvents > len(events) :
                print("WARNING: Trying to place",nEvents,"events from", fileName,
                      "but the file doesn't have enough. Rerolling, but this is not physical.")
                continue
            break
        if self.args.verbose : print ( "Placing",nEvents,"events from", fileName )

        # Get events 
        toPlace = self.rng.choice( a=events, size=nEvents, p=probs, replace=False )

        # Place at random times
        if not squash :
            times = list(self.rng.uniform(low=0, high=intTime, size=nEvents))        

        for inevt in toPlace :
            Time = 0
            if not squash :
                Time = times.pop()

            # print ( " -- Placing at", time)
            
            particles, vertices = [], []
            # Stores the vertices of the event inside a vertex container. These vertices are in increasing order so we can index them with [abs(vertex_id)-1]
            for vertex in inevt.vertices:
                position=vertex.position
                if not squash :
                    # Unit conversion
                    TimeHepmc = c_light*Time
                    position=position+pyhepmc.FourVector(x=0,y=0,z=0,t=TimeHepmc)
                v1=pyhepmc.GenVertex(position)
                vertices.append(v1)
            
            # copies the particles and attaches them to their corresponding vertices
            for particle in inevt.particles:
                # no copy/clone operator...
                momentum, status, pid = particle.momentum, particle.status, particle.pid                
                p1 = pyhepmc.GenParticle(momentum=momentum, pid=pid, status=status)
                p1.generated_mass = particle.generated_mass
                particles.append(p1)
                
                # since the beam particles do not have a production vertex they cannot be attached to a production vertex
                if particle.production_vertex.id < 0:
                    production_vertex=particle.production_vertex.id
                    vertices[abs(production_vertex)-1].add_particle_out(p1)
                    hepSlice.add_particle(p1)
                    
                # Adds particles with an end vertex to their end vertices
                if particle.end_vertex:
                    end_vertex = particle.end_vertex.id
                    vertices[abs(end_vertex)-1].add_particle_in(p1)

            # Adds the vertices with the attached particles to the event
            for vertex in vertices:
                hepSlice.add_vertex(vertex)

        return hepSlice

    # ============================================================================================
    def poissonTimes( self, mu, endTime ):
        """Return an np.array of poisson-distributed times."""
        #Exponential distribution describes the time between poisson. We could start with an array of expected length
        #and then cut or fill as needed. Not very readable for very little gain.
        t = 0
        ret=np.array([])
        while True :
            delt = self.rng.exponential( mu )
            t = t + delt
            if t >= endTime :
                break
            ret = np.append(ret,t)
        return ret
    
    # ============================================================================================
    def nameGen(self):
        # # datetime object containing current date and time
        # now = datetime.now()
        # dt_string = now.strftime("%Y_%m_%d_%H_%M_%S") # YY/mm/dd H:M:S        
        # return 'Sig_Back_Combo_{}_{}_event.hepmc'.format(dt_string,numEvents)
        name = self.args.signalFile
        if self.args.nSlices > 0 : 
            name = name.replace(".hepmc","_n_{}.hepmc".format(self.args.nSlices))
            
        name = "bgmerged_"+name
        return name        

    def squawk(self,i) :
        print('Working on slice {}'.format( i+1 ))
        print('Resident Memory',psutil.Process().memory_info().rss / 1024 / 1024,"MB")
    
    # ============================================================================================
    def fileWriter(combo_cont):
        numEvents=len(combo_cont)
        name=nameGen(numEvents)
        Event_ind=0
        with WriterAscii(name) as f:
            for i in combo_cont:
                i.event_number=Event_ind
                Event_ind=Event_ind+1
                f.write_event(i)
        f.close()

    
# ============================================================================================
# Main function
if __name__ == '__main__':

    t0 = time.time()
    sbm = signal_background_merger()
    sbm.merge()
    print()
    print('==================================================================')
    print('Overall running time:',np.round((time.time()-t0)/60.,2),'min')



