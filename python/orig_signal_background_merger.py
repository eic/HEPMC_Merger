import pyhepmc as hep
from pyhepmc.io import WriterAscii
import random as rd
from datetime import datetime
import argparse

# ============================================================================================
# Prints a banner
def banner(Signal_File,Background_Files,one_background_particle,int_window):
    print('\n==================================================================')
    print('*** EPIC HEPMC MERGER ***')
    print('author: Benjamen Sterwerf, UC Berkeley (bsterwerf@berkeley.edu)')
    print('\nRun this code as:\npython signal_background_merger.py --help\nto get more info')
    print('------------------------------------------------------------------')
    print('Signal events will be read from:')
    print('\t- ',Signal_File,'\n')
    print('Background files will be read from:')
    for bf in Background_Files:
        print('\t- ',bf)
    print('\none_background_particle:',one_background_particle)
    print('Int_Window:',int_window, 'ns')
    print('==================================================================\n')

# ============================================================================================
# opens HEPMC File in any generation and appends each event to a container which is then returned
def fileProcess(fileName):
    container=[]
    with hep.open(fileName) as f:
        for event in f:
            container.append(event)
            pass
    return container

# ============================================================================================
def merger(one_background_particle,Int_Window,sig_cont,back_cont):

    c_light  = 299.792458 # speed of light = 299.792458 mm/ns to get mm
    combo_cont=[] # container for combo events

    lengths = [len(bck) for bck in back_cont]
    lengths.extend([len(sig_cont)])

    # If files are provided with variable number of events, run as many events as events are in the smallest file
    # i here are events
    for i in range(min(lengths)):
        combo_cont.append(hep.GenEvent(hep.Units.GEV, hep.Units.MM))

        # -----------------------------------------------------
        # SIGNALS
        sig_particles, sig_vertices = [], []
        sig_time_shift = c_light*(Int_Window*rd.random())
        
        # Stores the vertices of the event inside a vertex container. These vertices are in increasing order so we can index them with [abs(vertex_id)-1]
        for vertex in sig_cont[i].vertices:
            position=vertex.position
            #if we are shifting the event in time uses the random time generator and adds this time to the position four vector for the vertex
            if Int_Window > 0:
                position=position+hep.FourVector(x=0,y=0,z=0,t=sig_time_shift)
            v1=hep.GenVertex(position)
            sig_vertices.append(v1)
        
        # copies the particles and attaches them to their corresponding vertices
        for particle in sig_cont[i].particles:
            momentum = particle.momentum
            status = particle.status
            pid = particle.pid
            
            # hep.GenParticle(momentum(px,py,pz,e), pdgid, status)
            p1 = hep.GenParticle(momentum, pid, status)
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
                position = position+hep.FourVector(x=0,y=0,z=0,t=back_time_shift)
                v1 = hep.GenVertex(position)
                back_vertices.append(v1)

            # copies the particles and attaches them to their corresponding vertices
            for particle in back_cont[b][i].particles:
                momentum = particle.momentum
                status = particle.status
                pid = particle.pid

                # hep.GenParticle(momentum(px,py,pz,e), pdgid, status)
                p1 = hep.GenParticle(momentum, pid, status)
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
def nameGen(numEvents):
    # datetime object containing current date and time
    now = datetime.now()
    dt_string = now.strftime("%Y_%m_%d_%H_%M_%S") # YY/mm/dd H:M:S
    return 'Sig_Back_Combo_{}_{}_event.hepmc'.format(dt_string,numEvents)

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

    parser = argparse.ArgumentParser(description='Run this code as: python signal_background_merger.py --Int_Window 2000 --one_background_particle False --Signal_File signal.hepmc --Background_Files background_1.hepmc background_2.hepmc')
    parser.add_argument('--one_background_particle', 
                        help='Set to True if all the background vertices in the file are to be connected temporally and spatially and False for backgrounds with many particles in the event like SR (synchrotron radiation)', 
                        action='store_true')
    parser.add_argument('--Int_Window', action='store', type=float, default=0.0,
                        help='length of the integration window in nanoseconds. Default is 0. If set to a positive value (in ns), the signal event will be moved to a random time within the integration window')
    parser.add_argument('--Signal_File', action='store',
                        help='Name of the HEPMC file with the signal events',default='dummy_signal.hepmc')
    parser.add_argument('--Background_Files', action='store', nargs='+',
                        help='Names of the HEPMC files with background events',default=['SR_out_single.hepmc.gz','dummy_background_2.hepmc'])
    
    args = parser.parse_args()    

    # Depending on the background type.
    # SR Synchrotron radiation will be set too false since there will be many individual photons
    one_background_particle = args.one_background_particle

    # opens up the signal file and extracts all individual events stores these events in sig_cont
    sig_cont=fileProcess(args.Signal_File)

    # opens up the background files and extracts all individual events. stores these events in back_cont
    back_cont = [ fileProcess(bck) for bck in args.Background_Files ]

    # integration window
    int_window = args.Int_Window

    banner(args.Signal_File,args.Background_Files,one_background_particle,int_window)

    # merger(one_background_particle,shift,Int_Window,sig_cont,back_cont)
    combo_cont= merger(args.one_background_particle,int_window,sig_cont,back_cont)

    fileWriter(combo_cont)
