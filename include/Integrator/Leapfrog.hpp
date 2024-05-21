#pragma once
#include "Integrator.hpp"

class Leapfrog : public Integrator {
    public:
    Leapfrog(double dt) : Integrator(dt) {}

    Particle integrate(const Particle& particle, double dt) override;
};