import numpy as np
import math
import random

class particle(object):
    '''
        :param position: A vector with the position of the particle in xyz coordinates in nm
        :param velocity: A vector with the velocity of the particle in xyz coordinates in nm/ps
        :param acceleration: A vector with the acceleration of the particle in xyz coordinates in nm/ps^2
        :param atom_number: Number of atoms in the model
        :param radius: Length of the radius of the bounding sphere in nm
        :param mass: Mass of the model in Daltons (this is equal to its molar mass)
        :param charge: Charge of the model in units of elementary charge
        :param pdb_id: The PDB id of the model
    '''
    def __init__(self, pdb_id="", atom_number=0, center=0.0, radius=0.0, mass=0.0, charge=0.0):
        self.position = np.array([0.0,0.0,0.0])
        self.velocity = np.array([0.0,0.0,0.0])
        self.acceleration = np.array([0.0,0.0,0.0])
        self.atom_number = atom_number
        self.center = center
        self.radius = radius
        self.mass = mass
        self.charge = charge
        self.pdb_id = pdb_id

class simulator(object):
    '''
        :param timestep: The time that is simulated between steps
        :param box_dimensions: A vector [xlength, ylength, zlength] with the dimensions of the box in nm. The bounds
                               for the coordinates are -length/2 < coordinate <= length/2, so that (0,0,0) is at the
                               center of the simulation box
        :param temperature: The temperature of the system in Kelvin
    '''
    def __init__(self, timestep, box_dimensions, temperature):
        self.timestep = timestep
        self.dimensions = box_dimensions
        self.temperature = temperature
        self.particles = []

    '''
        Returns the distance vector pointing from particle1 to particle2 according to periodic boundary conditions.
        For each dimension, if the distance in that dimension is more than half of the length of that dimension, we
        adjust the distance so that it is less than that.
    '''
    def get_distance(self, particle1, particle2):
        dist_vector = particle2.position - particle1.position
        
        # Adjust the x distance
        if dist_vector[0] > self.dimensions[0]/2: dist_vector[0] -= self.dimensions[0]
        if dist_vector[0] <= -self.dimensions[0]/2: dist_vector[0] += self.dimensions[0]

        # Adjust the y distance
        if dist_vector[1] > self.dimensions[1]/2: dist_vector[1] -= self.dimensions[1]
        if dist_vector[1] <= -self.dimensions[1]/2: dist_vector[1] += self.dimensions[1]
        
        # Adjust the z distance
        if dist_vector[2] > self.dimensions[2]/2: dist_vector[2] -= self.dimensions[2]
        if dist_vector[2] <= -self.dimensions[2]/2: dist_vector[2] += self.dimensions[2]
        
        return dist_vector

    '''
        Attempts to add a particle to the box. It chooses its xyz coordinates uniformly at random. It will keep choosing
        coordinates if there are collisions between other particles.
        
        :param new_particle: The particle to be added to the simulation
        :param num_attempts: The maximum number of attempts to insert the particle into the box. If it goes over this number,
            then it will not insert into the box. If this is -1, then it will never stop.
        
        :return: 0 if the particle was successfully inserted, -1 otherwise
    '''
    def add_particle(self, new_particle, num_attempts=-1):
        successful = False
        count = 0
        while not successful:
            # Check if this too many attempts
            count += 1
            if num_attempts > -1 and count > num_attempts:
                break

            # Set new coordinates for the particle
            new_particle.position[0] = random.uniform(-self.dimensions[0]/2, self.dimensions[0]/2)
            new_particle.position[1] = random.uniform(-self.dimensions[1]/2, self.dimensions[1]/2)
            new_particle.position[2] = random.uniform(-self.dimensions[2]/2, self.dimensions[2]/2)

            # Check for collisions
            successful = True
            for p in self.particles:
                dist = self.get_distance(new_particle, p)
                if np.sqrt(np.sum(np.power(dist, 2))) < new_particle.radius + p.radius:
                    successful = False
                    break

        if successful:
            self.particles.append(new_particle)
            return 0
        else:
            return -1

    '''
        Returns the force acted on particle1 by particle2 according to a Lennard-Jones potential

        :param energy_depth: The difference between the minimum potential and zero potential, in units of
                            nN * nm
        :param characteristic_radius: The distance at which the potential is zero, in units of nm

        :output: The force vector acted on particle1 by particle2 in units of nN
    '''
    def lj_force(self, particle1, particle2, energy_depth, characteristic_radius):
        dist_vector = self.get_distance(particle1, particle2)
        distance_sq = np.sum(np.power(dist_vector,2))
        distance = np.sqrt(distance_sq)
        term1 = particle1.radius*particle2.radius
        term2 = distance_sq - (particle1.radius + particle2.radius)**2
        term3 = distance_sq - (particle1.radius - particle2.radius)**2

        # Compute the magnitude of the attractive force
        attractive_force = -2*term1*(1/(term2**2) + 1/(term3**2))
        attractive_force += (1/term2 - 1/term3)
        attractive_force *= -4*energy_depth*distance*math.pi*math.pi/3

        # Compute the magnitude of the repulsive force
        repulsive_force = -(characteristic_radius**6)/((distance - particle1.radius - particle2.radius)**8)
        repulsive_force *= term1 / (particle1.radius + particle2.radius)
        repulsive_force *= 7*energy_depth*math.pi*math.pi/315

        # Return the vector pointing from particle1 to particle2 with magnitude equal to the sum of the two forces
        return ((attractive_force + repulsive_force) / distance) * dist_vector

    '''
        Returns the force acted on particle1 by particle2 according to a DLVO electrostatic potential

        :param debye_screening_length: The debye screening length in nm
        :param dielectric_constant: The relative permitivity of the medium

        :output: The force vector acted on particle1 by particle2 in units of nN
    '''
    def dlvo_electrostatic_force(self, particle1, particle2, dielectric_constant, debye_screening_length):
        kappa = 1/debye_screening_length
        dlvo_constant = 0.2307077556 # Equal to e^2/4*pi*vacuum permitivity, in units of nN * nm^2
        
        dist_vector = self.get_distance(particle1, particle2)
        distance = np.sqrt(np.sum(np.power(dist_vector, 2)))

        # Compute the force
        force = -math.exp(-kappa*(distance - particle1.radius - particle2.radius))/distance
        force *= 1/distance + kappa
        force *= particle1.charge*particle2.charge*dlvo_constant
        force /= dielectric_constant*(1 + kappa*particle1.radius)*(1 + kappa*particle2.radius)

        # Return the vector pointing from particle1 to particle2 with magnitude equal to the force
        return (force / distance) * dist_vector

    '''
        Computes the net force on each particle in the simulator's particle list.
        Force is in units of nN
    '''
    def compute_net_force(self):
        net_force = [np.array([0,0,0]) for x in range(len(self.particles))]
        for i in range(len(self.particles)):
            for j in range(i+1, len(self.particles)):
                # Compute the pairwise force between particles i and j
                force = self.lj_force(self.particles[i], self.particles[j], 1, 1)
                force += self.dlvo_electrostatic_force(self.particles[i], self.particles[j], 1, 1)
                
                # Add the force to the net force. We add the same force but in the opposite direction to j
                # according to Newton's 3rd law
                net_force[i] += force
                net_force[j] += -1*force
        
        return net_force

    '''
        Computes the new position of each particle. Should be the first method called in an update step.
        Positions are in nm.
    '''
    def update_positions(self):
        for i in range(len(self.particles)):
            self.particles[i].position += self.timestep*self.particles[i].velocity 
            self.particles[i].position += self.timestep * self.timestep * self.particles[i].acceleration / 2

            # Now make sure that the particle is within bounds
            if self.particles[i].position[0] > self.dimensions[0]/2: self.particles[i].position[0] -= self.dimensions[0]
            if self.particles[i].position[0] <= -self.dimensions[0]/2: self.particles[i].position[0] += self.dimensions[0]

            if self.particles[i].position[1] > self.dimensions[1]/2: self.particles[i].position[1] -= self.dimensions[1]
            if self.particles[i].position[1] <= -self.dimensions[1]/2: self.particles[i].position[1] += self.dimensions[1]

            if self.particles[i].position[2] > self.dimensions[2]/2: self.particles[i].position[2] -= self.dimensions[2]
            if self.particles[i].position[2] <= -self.dimensions[2]/2: self.particles[i].position[2] += self.dimensions[2]
        
    '''
        Computes the new acceleration of each particle. Should be called after updating positions.
        This does not set the particles acceleration to the new accelerations. 
        The accelerations are in nm/ps^2
    '''
    def get_new_acceleration(self):
        # This is avogadro's constant times 10^-24, and this is to convert nN to Da*nm/ps^2
        conversion_factor = 602.214076 
        net_force = self.compute_net_force()
        for i in range(len(self.particles)):
            net_force[i] = net_force[i] * conversion_factor / self.particles[i].mass
        
        return net_force

    '''
        Computes the new velocity of each particle. Should be called after getting the new accelerations.
        The velocities are in nm/ps. The ith new acceleration should correspond to particle i
    '''
    def update_velocities(self, new_accelerations):
        for i in range(len(self.particles)):
            self.particles[i].velocity += 0.5 * self.timestep * (self.particles[i].acceleration + new_accelerations[i])

    '''
        Updates the accelerations for each particle. Should be called after updating velocities
    '''
    def update_acceleration(self, new_accelerations):
        for i in range(len(self.particles)):
            self.particles[i].acceleration = new_accelerations[i]

    '''
        Adjusts the velocities so that the temperature is constant
    '''
    def adjust_velocities(self):
        # This is Boltzmann constant in units of Da*nm^2/(ps^2 * K)
        kb = 0.00831446261815324

        # Compute the two times the total kinetic energy in units of Da*nm^2/ps^2
        tot_kinetic_x = 0.0
        tot_kinetic_y = 0.0
        tot_kinetic_z = 0.0
        for p in self.particles:
            tot_kinetic_x += p.mass * p.velocity[0] * p.velocity[0]
            tot_kinetic_y += p.mass * p.velocity[1] * p.velocity[1]
            tot_kinetic_z += p.mass * p.velocity[2] * p.velocity[2]

        # Compute the velocity scaling factor in each dimension
        target_energy = len(self.particles) * kb * self.temperature
        scale_factor_x = math.sqrt(target_energy/tot_kinetic_x)
        scale_factor_y = math.sqrt(target_energy/tot_kinetic_y)
        scale_factor_z = math.sqrt(target_energy/tot_kinetic_z)

        # Update the velocities
        for p in self.particles:
            p.velocity[0] *= scale_factor_x
            p.velocity[1] *= scale_factor_y
            p.velocity[2] *= scale_factor_z

    '''
        Initalizes the velocities of the particles and also computes the accelerations of the particels
        Sets the total linear momentum to zero. At the moment, the velocities are assigned uniformly at
        random on the interval [-1, 1] before they are scaled

        TODO: Sample velocities from maxwell boltzmann distribution instead
    '''
    def init_simulation(self):
        avg_mom_x = 0.0
        avg_mom_y = 0.0
        avg_mom_z = 0.0

        # Assign the velocities and compute the average linear momentum
        for p in self.particles:
            p.velocity[0] = random.uniform(-1.0, 1.0)
            p.velocity[1] = random.uniform(-1.0, 1.0)
            p.velocity[2] = random.uniform(-1.0, 1.0)
            avg_mom_x += p.mass * p.velocity[0]
            avg_mom_y += p.mass * p.velocity[1]
            avg_mom_z += p.mass * p.velocity[2]

        avg_mom_x /= len(self.particles)
        avg_mom_y /= len(self.particles)
        avg_mom_z /= len(self.particles)

        # Adjust the velocities so that there is no net linear momentum
        for p in self.particles:
            p.velocity[0] -= avg_mom_x / p.mass
            p.velocity[1] -= avg_mom_y / p.mass
            p.velocity[2] -= avg_mom_z / p.mass

        # Scale the velocities so that they are at the desired temperature
        self.adjust_velocities()

        # Now assign the accelerations
        self.update_acceleration(self.get_new_acceleration())
    
    '''
        Performs one step of the simulation
    '''
    def single_step(self):
        # Update the position, velocity and acceleration
        self.update_positions()
        new_accelerations = self.get_new_acceleration()
        self.update_velocities(new_accelerations)
        self.update_acceleration(new_accelerations)

        # Adjust the velocities
        self.adjust_velocities()