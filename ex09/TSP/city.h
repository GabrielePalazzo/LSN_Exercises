#ifndef __city_h__
#define __city_h__

#define _CIRCLE_
#define CHINESE 0

#ifdef _BOX_
    #define numOfCities 32
    #define numOfPopulations 5000
    #define numOfGenerations 100
    //#define absoluteBest 4.63666
    #define absoluteBest 9.27331
    #define crossoverRate 0.7
    #define mutationRate 0.3
    #define minMutationRate 0.01
    #define maxMutationRate 0.09
#endif

#ifdef _CIRCLE_
    #define numOfCities 32
    #define numOfPopulations 5000
    #define numOfGenerations 100
    #define absoluteBest 6.2403
    #define crossoverRate 0.7
    #define mutationRate 0.3
    #define minMutationRate 0.1
    #define maxMutationRate 0.9
#endif

#ifdef _FILE_
    #define numOfCities 48
    #define numOfPopulations 16000
    #define numOfGenerations 100
    #define absoluteBest 12377.8
    #define crossoverRate 0.6
    #define mutationRate 0.2
    #define minMutationRate 0.06
    #define maxMutationRate 0.30
#endif


// (0.8 for the box)


class City {
public:
    City();
    City(double a, double b);
    ~City();
    double x;
    double y;
};

//double Old_Distance(City*, int*);
double Distance(double*, int*);
double Distance_2(double*, int*);
void Initialize_Distance(double*, City*);
double Dist(double*, int, int);
double Computed_Distance(City, City);

#endif