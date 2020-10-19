#include <iostream>
#include <cmath>
#include "metropolis.h"

using namespace std;

WaveFunction::~WaveFunction() { }

double Y_1_0_0::Probability(double r, double theta, double phi) {
    return exp(-r) * exp(-r) / M_PI;
}

double Y_2_1_0::Probability(double r, double theta, double phi) {
    //double self = (r * cos(theta) / sqrt(M_PI/2.0) / 8.0);
    //cout << theta << " ";
    //cout << exp(-r) * r * r / M_PI / 32.0 * cos(theta) * cos(theta) << " ";
    return exp(-r) * r * r / M_PI / 32.0 * cos(theta) * cos(theta);
}

void Write(double* array, ofstream& write) {
    for (int i = 0; i < numberDimensions; i++) write << array[i] << " ";
    write << endl;
}

void Metropolis_Step(double* currentPosition, Random& generator, int sampling, WaveFunction* wF, int& accepted) {
    double newPosition[numberDimensions];
    for (int i = 0; i < numberDimensions; i++) {
        if (sampling == UNIFORM)
            newPosition[i] = currentPosition[i] + (generator.Rannyu() - 0.5) * sigma;
        else if (sampling == GAUSSIAN)
            newPosition[i] = currentPosition[i] + generator.Gauss(0, sigma);
    }
    double rNew = 0, rOld = 0;
    for (int i = 0; i < numberDimensions; i++) {
        rNew += newPosition[i] * newPosition[i];
        rOld += currentPosition[i] * currentPosition[i];
    }
    rNew = sqrt(rNew);
    rOld = sqrt(rOld);
    
    double thetaNew = 0, thetaOld = 0;
    if (rNew > 0) thetaNew = acos(newPosition[2] / rNew);
    else thetaNew = 0.0; // I don't want to be stuck in (0,0,0) forever
    if (rOld > 0) thetaOld = acos(currentPosition[2] / rOld);
    else thetaOld = 0.0; // I don't want to be stuck in (0,0,0) forever
    
    double phiNew = 0, phiOld = 0;
    phiNew = atan(newPosition[1]/newPosition[0]);
    phiOld = atan(currentPosition[1]/currentPosition[0]);
    
    //cout << wF->Probability(rNew, thetaNew, phiNew) << " " << wF->Probability(rOld, thetaOld, phiOld) << endl;
    
    if (generator.Rannyu() < (wF->Probability(rNew, thetaNew, phiNew) / wF->Probability(rOld, thetaOld, phiOld))) {
        for (int i = 0; i < numberDimensions; i++) {
            currentPosition[i] = newPosition[i];
        }
        accepted++;
    }
    
    /*
    rOld = 0;
    for (int i = 0; i < numberDimensions; i++) {
        rOld += currentPosition[i] * currentPosition[i];
    }
    rOld = sqrt(rOld);
    
    cout << rOld << endl;*/
    
    //Write(currentPosition, write);
}
