#include "Particle.hpp"

Particle& Particle::operator=(const Particle& other){
    mass = other.mass;
    position = other.position;
    velocity = other.velocity;
    acceleration = other.acceleration;
    return *this;
}

Eigen::Vector3d Particle::getForceInteraction(const Particle& p1, const Particle& p2, double softening){
    const Vector3d dPosition = p2.position - p1.position;
    return G * p1.mass * p2.mass * dPosition / std::pow(dPosition.dot(dPosition) + softening, 1.5);
}

double Particle::getKineticEnergy(const Particle& p){
    return 0.5 * p.mass * p.velocity.dot(p.velocity);
}

double Particle::getPotentialEnergy(const Particle& p1, const Particle& p2, double softening){
    const Vector3d dPosition = p2.position - p1.position;
    return G * p1.mass * p2.mass / std::sqrt(dPosition.dot(dPosition) + softening);
}