import pyhepmc as hep
from pyhepmc.io import WriterAscii
import random as rd
from datetime import datetime
import argparse

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
def merger(one_background_particle,shift,Int_Window,sig_cont,back_cont):
    # container for combo events
    combo_cont=[]

    lengths = [len(bck) for bck in back_cont]
    lengths.extend([len(sig_cont)])

    # If files are provided with variable number of events, run as many events as events are in the smallest file
    # i here are events
    for i in range(min(lengths)):
        combo_cont.append(hep.GenEvent(hep.Units.GEV, hep.Units.MM))

        # -----------------------------------------------------
        # SIGNALS
        sig_particles=[]
        sig_vertices=[]
        sig_time_shift=c_light*(Int_Window*rd.random())
        # sig_vertices.append(hep.GenVertex(sig_cont[i].event_pos()))
        
        # Stores the vertices of the event inside a vertex container. These vertices are in increasing order so we can index them with [abs(vertex_id)-1]
        for vertex in sig_cont[i].vertices:
            position=vertex.position
            #if we are shifting the event in time uses the random time generator and adds this time to the position four vector for the vertex
            if shift:
                position=position+hep.FourVector(x=0,y=0,z=0,t=sig_time_shift)
            v1=hep.GenVertex(position)
            sig_vertices.append(v1)
        
        # copies the particles and attaches them to their corresponding vertices
        for particle in sig_cont[i].particles:
            momentum=particle.momentum
            status=particle.status
            pid=particle.pid
            
            # hep.GenParticle(momentum(px,py,pz,e), pdgid, status)
            p1 = hep.GenParticle(momentum, pid, status)
            p1.generated_mass = particle.generated_mass
            sig_particles.append(p1)
            # since the beam particles do not have a production vertex they cannot be attached to a production vertex
            if particle.production_vertex.id<0:
                production_vertex=particle.production_vertex.id
                sig_vertices[abs(production_vertex)-1].add_particle_out(p1)
                combo_cont[i].add_particle(p1)
            # Adds particles with an end vertex to their end vertices
            if particle.end_vertex:
                end_vertex=particle.end_vertex.id
                sig_vertices[abs(end_vertex)-1].add_particle_in(p1)
        # Adds the vertices with the attached particles to the event
        for vertex in sig_vertices:
            combo_cont[i].add_vertex(vertex)
        
        # -----------------------------------------------------
        # BACKGROUNDS
        # Loop over different background input files
        for b in range(len(back_cont)):

            back_particles=[]
            back_vertices=[]

            # Standard use except for SR backgrounds. Will only shift the background event when this conditional is triggered.
            if one_background_particle:
                #time to shift each background vertex
                back_time_shift=c_light*(Int_Window*rd.random())

            for vertex in back_cont[b][i].vertices:

                position=vertex.position
                # When the background event has many particles separated in time, like SR backgrounds, trigger this coniditional by changing the one_background_particle to False
                if one_background_particle==False:
                    back_time_shift=c_light*(Int_Window*rd.random())
                position=position+hep.FourVector(x=0,y=0,z=0,t=back_time_shift)
                v1=hep.GenVertex(position)
                back_vertices.append(v1)

            # copies the particles and attaches them to their corresponding vertices
            for particle in back_cont[b][i].particles:
                momentum=particle.momentum
                status=particle.status
                pid=particle.pid

                # hep.GenParticle(momentum(px,py,pz,e), pdgid, status)
                p1 = hep.GenParticle(momentum, pid, status)
                p1.generated_mass = particle.generated_mass
                back_particles.append(p1)
                # since the beam particles do not have a production vertex they cannot be attached to a production vertex
                if particle.production_vertex.id<0:
                    production_vertex=particle.production_vertex.id
                    back_vertices[abs(production_vertex)-1].add_particle_out(p1)
                    combo_cont[i].add_particle(p1)
                # Adds particles with an end vertex to their end vertices
                if particle.end_vertex:     
                    end_vertex=particle.end_vertex.id
                    back_vertices[abs(end_vertex)-1].add_particle_in(p1)
            # Adds the vertices with the attached particles to the event        
            for vertex in back_vertices:
                combo_cont[i].add_vertex(vertex)
    return combo_cont

# ============================================================================================
def nameGen(numEvents):
    # datetime object containing current date and time
    now = datetime.now()

    # YY/mm/dd H:M:S
    dt_string = now.strftime("%Y_%m_%d_%H_%M_%S")

    name='Sig_Back_Combo_{}_event_{}.hepmc'.format(numEvents,dt_string)
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

    parser = argparse.ArgumentParser(description='Merging Background and signal files')
    parser.add_argument('--one_background_particle', 
                        help='Set to True if all the background vertices in the file are to be connected temporally and spatially and False for backgrounds with many particles in the event like SR (synchrotron radiation)', 
                        action='store_true')
    parser.add_argument('--shifter', 
                        help='depending on how integration frames are decided, the signal event will either be positioned at time 0 (False) or set to a random time within the integration window', 
                        action='store_true')
    parser.add_argument('--Int_Window', action='store', type=int, default=2000,
                        help='length of the integration window in nanoseconds')
    parser.add_argument('--Signal_File', action='store',
                        help='Name of the HEPMC file with the signal events',default='Test_DIS_event.hepmc')
    parser.add_argument('--Background_Files', action='store', nargs='+',
                        help='Name of the HEPMC file with the background events',default='Test_Back_event.hepmc')
    
    #parser.add_argument('-l','--list', nargs='+', help='<Required> Set flag', required=True)

    args = parser.parse_args()

    print(args.Background_Files)

    # initiation of the variables
    c_light  = 299.792458 # speed of light = 299.792458 mm/ns to get mm

    # Whether or not we move the signal event
    shift=args.shifter
    shift=False

    # Depending on the background type.
    # SR Synchrotron radiation will be set too false since there will be many individual photons
    one_background_particle=args.one_background_particle

    print(one_background_particle)

    # opens up the signal file and extracts all individual events stores these events in sig_cont
    # with hep.open("pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1_000.hepmc") as f:
    sig_cont=fileProcess(args.Signal_File)

    # opens up the signal file and extracts all individual events. stores these events in back_cont
    # with hep.open("photon_event16_03_2023_02_05_02.hepmc") as f:
    #back_cont=fileProcess(args.Background_Files)
    back_cont = [ fileProcess(bck) for bck in args.Background_Files ]

    # merger(one_background_particle,shift,Int_Window,sig_cont,back_cont)
    combo_cont= merger(args.one_background_particle,args.shifter,args.Int_Window,sig_cont,back_cont)

    fileWriter(combo_cont)