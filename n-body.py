import numpy as np
import matplotlib.pyplot as plt

class Particle:
	""" Class representing a particle in a 2D space
        
    :param x: x-coordinate of the particle
	:param y: y-coordinate of the particle
	:param vx: x-component of the velocity of the particle
	:param vy: y-component of the velocity of the particle
	:param mass: mass of the particle
    """

	def __init__(self, x: float, y: float, vx: float, vy: float, mass: float):
		self.position = np.array([x, y])
		self.velocity = np.array([vx, vy])
		self.mass = mass
		
	def getAcceleration(self, rhs: 'Particle', softening: float) -> np.array:
		""" Returns the acceleration of the particle due to another particle
		
		:param rhs: the other particle
		:param softening: softening factor to avoid division by zero
		"""
		G = 6.67430e-11
		dPosition = rhs.position - self.position
		return G * rhs.mass * dPosition / (np.dot(dPosition, dPosition) + softening)**1.5

	def getPotentialEnergy(self, rhs: 'Particle', softening: float) -> float:
		""" Returns the potential energy of the particle due to another particle

		:param rhs: the other particle
		:param softening: softening factor to avoid division by zero
		"""
		G = 6.67430e-11
		dPosition = rhs.position - self.position
		return G * self.mass * rhs.mass / np.sqrt(np.dot(dPosition, dPosition) + softening)
	
	def getKineticEnergy(self) -> float:
		""" Returns the kinetic energy of the particle """
		return 0.5 * self.mass * np.dot(self.velocity, self.velocity)

	def move(self, acceleration: np.array, dt: float) -> None:
		""" Moves the particle based on the acceleration and the time step
		:param acceleration: the acceleration acting on the particle
		:param dt: the time step
		"""
		self.velocity += acceleration * dt
		self.position += self.velocity * dt
		
class StarSystem:
	""" Class representing a system of particles. The softening parameter is calculated based on this paper: https://academic.oup.com/mnras/article/314/3/475/969154

	:param particles: list of particles in the system
	"""

	def __init__(self, particles: list[Particle]):
		self.particles = particles
		self.n = len(particles)
		self.softening = (0.98 * self.n**-0.26)**2
		
	def getAccelerations(self):
		""" Returns the accelerations of all particles in the system """
		accelerations = []
		for i in range(self.n):
			acceleration = np.zeros(2)
			for j in range(self.n):
				if i != j:
					acceleration += self.particles[i].getAcceleration(self.particles[j], self.softening)
			accelerations.append(acceleration)

		return accelerations

	def kineticEnergy(self):
		""" Returns the total kinetic energy of the system """
		energy = 0
		for particle in range(self.particles):
			energy += particle.getKineticEnergy()
		return energy

	def potentialEnergy(self):
		""" Returns the total potential energy of the system """
		energy = 0
		for i in range(self.n):
			for j in range(i+1, self.n):
				energy += self.particles[i].getPotentialEnergy(self.particles[j], self.softening)
		return energy
	
	def getTotalEnergy(self):
		""" Returns the total energy of the system """
		return self.kineticEnergy() + self.potentialEnergy()

	def evolve(self, dt: float):
		""" Evolves the system based on the time step """
		accelerations = self.getAccelerations()
		for i in range(self.n):
			self.particles[i].move(accelerations[i], dt)
	
def solarSystem():
	""" Returns a list of particles representing the solar system """
	sun     = Particle(0.0,      0.0,      0.0,    0.0,      1.98847e30)
	mercury = Particle(5.79e10,  0.0,      0.0,    4.74e4,   3.3011e23)
	venus   = Particle(0.0,      1.08e11,  3.50e4, 0,        4.8673e24)
	earth   = Particle(-1.50e11, 0,        0.0,    2.9783e4, 5.97219e24)
	mars    = Particle(0.0,      -2.28e11, 2.41e4, 0,        6.4169e23)
	jupiter = Particle(7.78e11,  0,        0.0,    1.31e4,   1.8981e27)
	saturn  = Particle(0.0,      1.43e12,  9.64e3, 0,        5.68232e26)
	uranus  = Particle(-2.88e12, 0.0,      0.0,    6.7991e3, 8.6810e25)
	neptune = Particle(0.0,      -4.50e12, 5.43e3, 0,        1.0241e26)
	return [sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]

def main():
	DAY = 24*60*60
	MONTH = 30*DAY
	YEAR = 12*MONTH

	tEnd = 100 * YEAR  
	dt = 5 * DAY
	iterations = int(tEnd/dt) 

	system = StarSystem(solarSystem())

	ax1 = plt.subplot(111)	
	ax1.set_title('Solar System')
	ax1.set_xlabel('x (m)')
	ax1.set_ylabel('y (m)')
	
	mercuryTrajectory = ax1.plot([], [], color='gray')
	venusTrajectory = ax1.plot([], [], color='orange')
	earthTrajectory = ax1.plot([], [], color='blue')
	marsTrajectory = ax1.plot([], [], color='red')
	jupiterTrajectory = ax1.plot([], [], color='brown')
	saturnTrajectory = ax1.plot([], [], color='brown')
	uranusTrajectory = ax1.plot([], [], color='blue')
	neptuneTrajectory = ax1.plot([], [], color='blue')
	
	plt.sca(ax1)
	ax1.set(xlim=(-5e12, 5e12), ylim=(-5e12, 5e12))
	ax1.set_aspect('equal', 'box')

	for _ in range(iterations):
		system.evolve(dt)
		
		mercuryTrajectory[0].set_data(
			np.append(mercuryTrajectory[0].get_xdata(), system.particles[1].position[0]),
			np.append(mercuryTrajectory[0].get_ydata(), system.particles[1].position[1]))

		venusTrajectory[0].set_data(
			np.append(venusTrajectory[0].get_xdata(), system.particles[2].position[0]),
			np.append(venusTrajectory[0].get_ydata(), system.particles[2].position[1]))

		earthTrajectory[0].set_data(
			np.append(earthTrajectory[0].get_xdata(), system.particles[3].position[0]),
			np.append(earthTrajectory[0].get_ydata(), system.particles[3].position[1]))

		marsTrajectory[0].set_data(
			np.append(marsTrajectory[0].get_xdata(), system.particles[4].position[0]),
			np.append(marsTrajectory[0].get_ydata(), system.particles[4].position[1]))

		jupiterTrajectory[0].set_data(
			np.append(jupiterTrajectory[0].get_xdata(), system.particles[5].position[0]),
			np.append(jupiterTrajectory[0].get_ydata(), system.particles[5].position[1]))

		saturnTrajectory[0].set_data(
			np.append(saturnTrajectory[0].get_xdata(), system.particles[6].position[0]),
			np.append(saturnTrajectory[0].get_ydata(), system.particles[6].position[1]))

		uranusTrajectory[0].set_data(
			np.append(uranusTrajectory[0].get_xdata(), system.particles[7].position[0]),
			np.append(uranusTrajectory[0].get_ydata(), system.particles[7].position[1]))

		neptuneTrajectory[0].set_data(
			np.append(neptuneTrajectory[0].get_xdata(), system.particles[8].position[0]),
			np.append(neptuneTrajectory[0].get_ydata(), system.particles[8].position[1]))

		plt.pause(0.000001)	
		plt.draw()
  
if __name__== "__main__":
  main()