#ifndef __VMC_H__
#define __VMC_H__

#include "metropolis.h"

void Minimize(WaveFunction*, double, Random&, double*, int, double*, double&);
void Simulate_Annealing(WaveFunction*, double, Random&, double*, int, double*, double&, double);
void Gradient_Descent(WaveFunction*, double, Random&, double*, int, double*, double&);

#endif
