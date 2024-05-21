#pragma once
#include "Particle/Particle.hpp"

class Integrator {
    public:
    Integrator(double dt) : dt(dt) {}

    virtual ~Integrator() = default;
    
    virtual Particle integrate(const Particle& particle, double dt) = 0;

    protected:
    double dt;
};