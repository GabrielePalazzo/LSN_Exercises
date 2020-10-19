#ifndef __metropolis_h__
#define __metropolis_h__

#include <fstream>
#include "random.h"

using namespace std;

#define numberDimensions 1
#define UNIFORM 0
#define GAUSSIAN 1

#define numberOfSteps 100000
#define numberOfBlocks 200
#define samplingType UNIFORM
#define equilibriumSteps 10000
#define MEASURE 1
#define MINIMIZE 0

//#define MU -0.78268
//#define SIGMA 0.659058
#define MU -0.78268
#define SIGMA 0.659058
#define delta 5.0

class WaveFunction {
public:
    //WaveFunction() { }
    virtual ~WaveFunction() = 0;
    virtual double Probability(double) = 0;
    virtual double Energy(double) = 0;
    virtual double Grad_sigma(double) = 0;
    virtual double Grad_mu(double) = 0;
    virtual void Set_mu(double) = 0;
    virtual void Set_sigma(double) = 0;
    virtual double Get_mu() = 0;
    virtual double Get_sigma() = 0;
};

class Variational: public WaveFunction {
private:
    double _mu;
    double _sigma;
public:
    Variational(double, double);
    ~Variational() { }
    double Probability(double);
    double Energy(double);
    double Grad_sigma(double);
    double Grad_mu(double);
    void Set_mu(double a) {_mu = a;}
    void Set_sigma(double a) {_sigma = a;}
    double Get_mu() {return _mu;}
    double Get_sigma() {return _sigma;}
};

void Write(double, ofstream&);
void Write(double, double, ofstream&);
void Metropolis_Step(double&, Random&, int, WaveFunction*, int&);

#endif
