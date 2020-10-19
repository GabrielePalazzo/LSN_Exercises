#include <iostream>
#include <cmath>
#include "metropolis.h"

using namespace std;

WaveFunction::~WaveFunction() { }

Variational::Variational(double a, double b) {
    _mu = a;
    _sigma = b;
}

double Variational::Probability(double r) {
    double denominator = 2 * _sigma * sqrt(M_PI) * (1 + exp(-(_mu / _sigma) * (_mu / _sigma)));
    return (exp(-(r - _mu) / _sigma * (r - _mu) / _sigma)
        + exp(-(r + _mu) / _sigma * (r + _mu) / _sigma)
        + 2* exp(-(r * r + _mu * _mu) / _sigma / _sigma)) / denominator;
}

double Variational::Energy(double x) {
    return (_sigma * _sigma - x * x - _mu * _mu + 2 * x * _mu * tanh(x * _mu / (_sigma * _sigma)))
        / (2 * _sigma * _sigma * _sigma * _sigma) + x * x * x * x - 2.5 * x * x;
}

double Variational::Grad_mu(double x) {
    return (- _mu + x * tanh(x * _mu / (_sigma * _sigma)) + x * _mu * (1 - tanh(x * _mu / (_sigma * _sigma)) * tanh(x * _mu / (_sigma * _sigma))) * x / (_sigma * _sigma))
    / (_sigma * _sigma * _sigma * _sigma);
}

double Variational::Grad_sigma(double x) {
    return (-2.0 / _sigma * (_sigma * _sigma - x * x - _mu * _mu + 2 * x * _mu * tanh(x * _mu / (_sigma * _sigma)))
            + (_sigma - 2 * x * _mu * (1 - tanh(x * _mu / (_sigma * _sigma)) * tanh(x * _mu / (_sigma * _sigma))) * x * _mu / (_sigma * _sigma * _sigma)))
    / (_sigma * _sigma * _sigma * _sigma);
    
    //return (_sigma * _sigma - x * x - _mu * _mu + 2 * x * _mu * tanh(x * _mu / (_sigma * _sigma)))
    // / (2 * _sigma * _sigma * _sigma * _sigma) + x * x * x * x - 2.5 * x * x;
}

void Write(double value, ofstream& write) {
    write << value << endl;
}

void Write(double value, double anotherValue, ofstream& write) {
    write << value << " " << anotherValue << endl;
}

void Metropolis_Step(double& currentPosition, Random& generator, int sampling, WaveFunction* wF, int& accepted) {
    double newPosition = 0.0;
    if (sampling == UNIFORM)
        newPosition = currentPosition + (generator.Rannyu() - 0.5) * delta;
    else if (sampling == GAUSSIAN)
        newPosition = currentPosition + generator.Gauss(0, delta);
    

    //cout << wF->Probability(rNew, mu, sigma) << " " << wF->Probability(rOld, mu, sigma) << endl;
    if (generator.Rannyu() < (wF->Probability(newPosition) / wF->Probability(currentPosition))) {
        currentPosition = newPosition;
        accepted++;
    }
    //cout << currentPosition[0] << endl;
    //cout << wF->Energy(currentPosition[0], mu, sigma) << endl;
    
    /*
    rOld = 0;
    for (int i = 0; i < numberDimensions; i++) {
        rOld += currentPosition[i] * currentPosition[i];
    }
    rOld = sqrt(rOld);
    
    cout << rOld << endl;*/
    
    //Write(currentPosition, write);
}
