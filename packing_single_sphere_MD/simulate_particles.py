import numpy as np
import math
'''
class particle(object):
    def __init__(self, x, y, z, atom_number, radius, mass, charge, pdb_id):
        self.position = np.array([x,y,z]) # Global coordinates of the center of the sphere
        self.atom_number = atom_number # Number of atoms in the model
        self.radius = radius # Radius of the bounding sphere
        self.mass = mass # Mass of the model
        self.charge = charge # Total charge of the model
        self.id = pdb_id # The pdb id of the model


class simulator(object):
    def __init__(self, timestep, number_of_steps):
        self.timestep = timestep
        self.number_of_steps = number_of_steps
        self.particles = []
        self.intermolecular_forces = []
'''

'''
    Particles are dictionaries in the format of
    {
        'pdb_id':      the pdb id of the particle
        'position':    numpy vector of x,y,z coordinates, in nm
        'atom_number': total number of atoms in the particle
        'radius':      radius of the bounding sphere of the particle, in nm
        'mass':        mass of the particle in Daltons (this is just its molar mass)
        'charge':      charge of the particle in units of elementary charge
    }
'''

'''
    Returns the force acted on particle1 by particle2 according to a Lennard-Jones potential

    :param energy_depth: The difference between the minimum potential and zero potential, in units of
                         nN * nm
    :param characteristic_radius: The distance at which the potential is zero, in units of nm

    :output: The force vector acted on particle1 by particle2 in units of nN
'''
def lj_force(particle1, particle2, energy_depth, characteristic_radius):
    dist_vector = particle1['position'] - particle2['position']
    distance_sq = np.sum(np.power(dist_vector,2))
    distance = np.sqrt(distance_sq)
    term1 = particle1['radius']*particle2['radius']
    term2 = distance_sq - (particle1['radius'] + particle2['radius'])**2
    term3 = distance_sq - (particle1['radius'] - particle2['radius'])**2

    # Compute the magnitude of the attractive force
    attractive_force = -2*term1*(1/(term2**2) + 1/(term3**2))
    attractive_force += (1/term2 - 1/term3)
    attractive_force *= -4*energy_depth*distance*math.pi*math.pi/3

    # Compute the magnitude of the repulsive force
    repulsive_force = -(characteristic_radius**6)/((distance - particle1['radius'] - particle2['radius'])**8)
    repulsive_force *= term1 / (particle1['radius'] + particle2['radius'])
    repulsive_force *= 7*energy_depth*math.pi*math.pi/315

    # Return the vector pointing from particle1 to particle2 with magnitude equal to the sum of the two forces
    return ((attractive_force + repulsive_force) / distance) * dist_vector

'''
    Returns the force acted on particle1 by particle2 according to a DLVO electrostatic potential

    :param debye_screening_length: The debye screening length in nm
    :param dielectric_constant: The relative permitivity of the medium

    :output: The force vector acted on particle1 by particle2 in units of nN
'''
def dlvo_electrostatic_force(particle1, particle2, dielectric_constant, debye_screening_length):
    kappa = 1/debye_screening_length
    dlvo_constant = 0.2307077556 # Equal to e^2/4*pi*vacuum permitivity, in units of nN * nm^2
    
    dist_vector = particle1['position'] - particle2['position']
    distance = np.sqrt(np.sum(np.power(dist_vector, 2)))

    # Compute the force
    force = -math.exp(-kappa*(distance - particle1['radius'] - particle2['radius']))/distance
    force *= 1/distance + kappa
    force *= particle1['charge']*particle2['charge']*dlvo_constant
    force /= dielectric_constant*(1 + kappa*particle1['radius'])*(1 + kappa*particle2['radius'])

    # Return the vector pointing from particle1 to particle2 with magnitude equal to the force
    return (force / distance) * dist_vector

'''
    Returns a vector of the updated position using Verlet integration

    :param prev_position: The position of the particle one timestep before its current position, in nm
    :param net_force: The vector of the net force acting on the particle, in nN
    :param timestep: The amount of time to increment, in ps

    :output: A vector representing the position the particle would be in if a constant net_force acted on 
             the particle for timestep unit of time
'''
def get_new_position(particle, prev_position, net_force, timestep):
    conversion_factor = 602.214076 # This is avogadro's constant times 10^-24, and this is to convert nN to Da*nm/ps^2
    acceleration = net_force * conversion_factor / particle['mass']
    
    return 2*particle['position'] - prev_position + acceleration * timestep * timestep

'''
    @param 
'''
def simulate_step():
    pass

'''
    @param N: the number of molecules to simulate
    @param 
'''
def simulate():
    pass