#include "Integrator/Leapfrog.hpp"

Particle Leapfrog::integrate(const Particle& particle, double dt){
    using Vector3d = Eigen::Vector3d;

    const Vector3d position = particle.getPosition() + 
                                particle.getVelocity() * dt + 
                                particle.getAcceleration() * 0.5 * dt * dt;
    const Vector3d velocity = particle.getVelocity() + particle.getAcceleration() * dt;

    return Particle(particle.getMass(), position, velocity, particle.getAcceleration());
}