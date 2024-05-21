#pragma once
#include "Integrator.hpp"

class Euler : public Integrator {
    public:
    Euler(double dt) : Integrator(dt) {}

    Particle integrate(const Particle& particle, double dt) override;
};