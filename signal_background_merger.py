import pyhepmc
from pyhepmc.io import WriterAscii
import random as rd
from datetime import datetime
import argparse

import numpy as np
import time
import sys

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
        
        parser.add_argument('-i','--signalFile', default='dummy_signal.hepmc',
                            help='Name of the HEPMC file with the signal events')
        parser.add_argument('-sf','--signalFreq', type=float, default=0.0,
                            help='Poisson-mu of the signal frequency in ns. Default is 0 to have exactly one signal event per slice. Set to the estimated DIS to randomize.')

    
        parser.add_argument('-b1','--bg1File', default='dummy_bg1.hepmc',
                            help='Name of the first HEPMC file with background events')
        parser.add_argument('-f1','--bg1Freq', type=float, default=1852.0,
                            help='Poisson-mu of the first background frequency in ns. Default is the estimated DIS frequency of 1852 ns. Set to 0 to use the weights in the corresponding input file.')

        parser.add_argument('-b2','--bg2File', default='dummy_bg2.hepmc',
                            help='Name of the second HEPMC file with background events')
        parser.add_argument('-f2','--bg2Freq', type=float, default=1852.0,
                            help='Poisson-mu of the second background frequency in ns. Default is the estimated DIS frequency of 1852 ns. Set to 0 to use the weights in the corresponding input file.')

        parser.add_argument('-b3','--bg3File', default='SR_out_single.hepmc',
                            help='Name of the third HEPMC file with background events')
        parser.add_argument('-f3','--bg3Freq', type=float, default=0,
                            help='Poisson-mu of the third background frequency in ns. Default is 0 to use the weights in the corresponding input file. Set to a value >0 to specify a poisson mu instead.')


        parser.add_argument('-o','--outputFile', default='',
                            help='Specify the output file name. By default it will be auto-generated.')

        parser.add_argument('-w','--intWindow', type=float, default=2000.0,
                            help='Length of the integration window in nanoseconds. Default is 2000.')
        parser.add_argument('-N','--nSlices', type=int, default=10,
                            help='Number of sampled time slices ("events"). Default is 10. If set to -1, all events in the signal file will be used and background files cycled as needed.')
        parser.add_argument('--squashTime', action='store_true',
                            help='Integration is performed but no time information is associated to vertices.')

        parser.add_argument('--rngSeed', action='store',type=int, default=None,
                            help='Random seed, default is None')

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

        # Open signal file
        try :
            sigFile=pyhepmc.io.ReaderAscii(self.args.signalFile) 
        except IOError as e:
            print ('Opening files failed: %s' % e.strerror)
            sys.exit()

        # Open input file readers, group them with the associated frequency,
        # then sort them into signal, frequency bg, or weight type bg
        self.freqBundles=[]
        self.weightBundles=[]
        self.makeBundles ( self.args.signalFile, self.args.signalFreq, signal=True )
        self.makeBundles ( self.args.bg1File, self.args.bg1Freq )
        self.makeBundles ( self.args.bg2File, self.args.bg2Freq )
        self.makeBundles ( self.args.bg3File, self.args.bg3Freq )
        
        # Open the output file
        if self.args.outputFile != "" :
            outputFileName = self.args.outputFile
        else :
            outputFileName = self.nameGen()

        # Process and write
        nSlices=self.args.nSlices
        with WriterAscii(outputFileName) as f:
            i=0
            while True :
                if (i==nSlices) :
                    print ( "Finished all requested slices." )
                    break
                if  i % 10000 == 0 : print('Working on slice {}'.format( i+1 ))
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


        # Clean up, close all input files
        
        sys.exit()
    
    # ============================================================================================
    def makeBundles( self, fileName, freq, signal=False ):
        """Create background timeline chunks, open input file, sort into frequency or weight type"""
        if fileName == "" : return

        try :
            # Should use File=pyhepmc.open(fileName) but that doesn't currently work
            File=pyhepmc.io.ReaderAscii(fileName)
        except IOError as e:
            print ('Opening {} failed: {}', fileName, e.strerror)
            sys.exit()

        if signal :
            bundle = [ File, freq ]
            self.sigBundle = bundle

        if freq<=0 :
            # file has its own weights, nothing else to do
            self.weightBundles.append(File)
            return

        bundle = [ File, freq ]
        self.freqBundles.append(bundle)

        return

    # ============================================================================================
    def mergeSlice(self, i):
        """Arrange the composition of an individual time slice"""
        
        hepSlice = pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM)
        
        # Signal and frequency background are handled very similarly,
        # handle them in one method and just ise a flag for what little is special
        try :
            hepSlice = self.addFreqEvents( self.sigBundle, hepSlice, signal=True )
        except EOFError :
            return
        
        # Treat frequency backgrounds very similarly 
        for freqBundle in self.freqBundles :
            freqFile, freqFreq = freqBundle
            hepSlice = self.addFreqEvents( freqBundle, hepSlice )

            

        return hepSlice
        


    # If files are provided with variable number of events, run as many events as events are in the smallest file
        # i here are events
        for i in range(min(lengths)):
            combo_cont.append(pyhepmc.GenEvent(pyhepmc.Units.GEV, pyhepmc.Units.MM))
            
            # -----------------------------------------------------
            # SIGNALS
            
            # Stores the vertices of the event inside a vertex container. These vertices are in increasing order so we can index them with [abs(vertex_id)-1]
            for vertex in sig_cont[i].vertices:
                position=vertex.position
                #if we are shifting the event in time uses the random time generator and adds this time to the position four vector for the vertex
                if Int_Window > 0:
                    position=position+pyhepmc.FourVector(x=0,y=0,z=0,t=sig_time_shift)
                v1=pyhepmc.GenVertex(position)
                sig_vertices.append(v1)
            
            # copies the particles and attaches them to their corresponding vertices
            for particle in sig_cont[i].particles:
                momentum = particle.momentum
                status = particle.status
                pid = particle.pid
                
                # pyhepmc.GenParticle(momentum(px,py,pz,e), pdgid, status)
                p1 = pyhepmc.GenParticle(momentum, pid, status)
                p1.generated_mass = particle.generated_mass
                sig_particles.append(p1)
                # since the beam particles do not have a production vertex they cannot be attached to a production vertex
                if particle.production_vertex.id < 0:
                    production_vertex=particle.production_vertex.id
                    sig_vertices[abs(production_vertex)-1].add_particle_out(p1)
                    combo_cont[i].add_particle(p1)
                # Adds particles with an end vertex to their end vertices
                if particle.end_vertex:
                    end_vertex = particle.end_vertex.id
                    sig_vertices[abs(end_vertex)-1].add_particle_in(p1)

            # Adds the vertices with the attached particles to the event
            for vertex in sig_vertices:
                combo_cont[i].add_vertex(vertex)
        
            # -----------------------------------------------------
            # BACKGROUNDS
            # Loop over different background input files
            for b in range(len(back_cont)):
                back_particles, back_vertices = [], []

                # Standard use except for SR backgrounds. Will only shift the background event when this conditional is triggered.
                if one_background_particle:
                    #time to shift each background vertex
                    back_time_shift = c_light*(Int_Window*rd.random())

                for vertex in back_cont[b][i].vertices:
                    position=vertex.position
                    # When the background event has many particles separated in time, like SR backgrounds, trigger this coniditional by changing the one_background_particle to False
                    if one_background_particle == False:
                        back_time_shift = c_light*(Int_Window*rd.random())
                    position = position+pyhepmc.FourVector(x=0,y=0,z=0,t=back_time_shift)
                    v1 = pyhepmc.GenVertex(position)
                    back_vertices.append(v1)
                    
                # copies the particles and attaches them to their corresponding vertices
                for particle in back_cont[b][i].particles:
                    momentum = particle.momentum
                    status = particle.status
                    pid = particle.pid

                    # pyhepmc.GenParticle(momentum(px,py,pz,e), pdgid, status)
                    p1 = pyhepmc.GenParticle(momentum, pid, status)
                    p1.generated_mass = particle.generated_mass
                    back_particles.append(p1)
                    # since the beam particles do not have a production vertex they cannot be attached to a production vertex
                    if particle.production_vertex.id < 0:
                        production_vertex = particle.production_vertex.id
                        back_vertices[abs(production_vertex)-1].add_particle_out(p1)
                        combo_cont[i].add_particle(p1)
                    # Adds particles with an end vertex to their end vertices
                    if particle.end_vertex:     
                        end_vertex = particle.end_vertex.id
                        back_vertices[abs(end_vertex)-1].add_particle_in(p1)
                # Adds the vertices with the attached particles to the event        
                for vertex in back_vertices:
                    combo_cont[i].add_vertex(vertex)
        return combo_cont

    # ============================================================================================
    def addFreqEvents( self, Bundle, hepSlice, signal=False ):
        """Handles the signal as well as frequency=style backgrounds"""
        
        File, Freq = Bundle
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
            slice = np.array([ self.rng.uniform(0,intTime) ])
        else:
            # Generate poisson-distributed times to place events
            slice = self.poissonTimes( Freq, intTime )

        print ( " !!! Adding {} events".format(slice.size) )
        if slice.size == 0 : return hepSlice

        # Insert events at all specified locations
        for Time in slice :
            ### Arrgh, GenEvent==None throws an exception
            try : 
                inevt = File.read()
                if inevt==None :
                    if signal :
                        # Exhausted signal events
                        raise EOFError
                    else :
                        # background file reached its end, reset to the start
                        ## TODO: do that...
                        print("hello bg error")
                        raise EOFError
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
    print('Overall running time:',np.round((time.time()-t0)/60.,2),'min')



