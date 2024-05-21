#pragma once
#include <Eigen/Core>

class Particle{
    using Vector3d = Eigen::Vector3d;

    public:
    Particle(double mass, 
             const Eigen::Vector3d& position, 
             const Eigen::Vector3d& velocity, 
             const Eigen::Vector3d& acceleration)
        : mass(mass), position(position), velocity(velocity), acceleration(acceleration) {}

    Particle& operator=(const Particle& other);

    inline double getMass() const{ 
        return mass; 
    }

    inline const Vector3d& getPosition() const{ 
        return position; 
    }

    inline const Vector3d& getVelocity() const{ 
        return velocity; 
    }

    inline const Vector3d& getAcceleration() const{ 
        return acceleration; 
    }

    static Vector3d getForceInteraction(const Particle& p1, const Particle& p2, double softening);

    static double getKineticEnergy(const Particle& p);

    static double getPotentialEnergy(const Particle& p1, const Particle& p2, double softening);

    private:
    double mass;
    Vector3d position;
    Vector3d velocity;
    Vector3d acceleration;

    const static double G = 6.67430e-11;
};