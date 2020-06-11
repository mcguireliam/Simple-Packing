import simulator
import pdb2ball_single
import random_select

'''
    Converts a protein dictionary to a particle
    
    :param pdbdict: A dictionary in the following format
        {'pdb_id': ****, 'atom_number': ****, 'center': ****, 'radius': ****, 'charge': ****, 'mass': ****}
'''
def pdbdict2particle(pdbdict):
    output = simulator.particle()
    output.atom_number = pdbdict['atom_number']
    output.center = pdbdict['center']
    output.radius = pdbdict['radius']
    output.charge = pdbdict['charge']
    output.mass = pdbdict['mass']
    output.pdb_id = pdbdict['pdb_id']

    return output

'''
    Attempts to add a particle that corresponds to protein to the simulation

    :param protein: A dictionary in the following format
        {'pdb_id': ****, 'atom_number': ****, 'center': ****, 'radius': ****, 'charge': ****, 'mass': ****}
    :param num_attempts: Number of attempts to insert the protein before not inserting it

    :return: 0 if insertion was successful, -1 otherwise
'''
def add_protein_to_simulation(sim: simulator.simulator, protein, num_attempts=-1):
    p = pdbdict2particle(protein)

    return sim.add_particle(p, num_attempts=num_attempts)