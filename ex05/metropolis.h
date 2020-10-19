#ifndef __metropolis_h__
#define __metropolis_h__

#include <fstream>
#include "random.h"

using namespace std;

#define numberDimensions 3
#define UNIFORM 0
#define GAUSSIAN 1
#define GROUND 0
#define EXCITED 1

#define numberOfSteps 1000000
#define samplingType UNIFORM
#define state EXCITED
#define equilibriumSteps 80000

#if samplingType == UNIFORM
    #if state == GROUND
        #define sigma 2.4
    #elif state == EXCITED
        #define sigma 5.8
    #else
        #define sigma 2.0
    #endif
#elif samplingType == GAUSSIAN
    #if state == GROUND
        #define sigma 0.75
    #elif state == EXCITED
        #define sigma 1.8
    #else
        #define sigma 0.8
    #endif
#else
    #define sigma 0.0
#endif

class WaveFunction {
public:
    //WaveFunction() { }
    virtual ~WaveFunction() = 0;
    virtual double Probability(double, double, double) = 0;
};

class Y_1_0_0: public WaveFunction {
public:
    Y_1_0_0() { }
    ~Y_1_0_0() { }
    double Probability(double, double, double);
};

class Y_2_1_0: public WaveFunction {
public:
    Y_2_1_0() { }
    ~Y_2_1_0() { }
    double Probability(double, double, double);
};

void Write(double*, ofstream&);
void Metropolis_Step(double*, Random&, int, WaveFunction*, int&);

#endif
