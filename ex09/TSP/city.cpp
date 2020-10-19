#include <cmath>
#include "city.h"

City::City () {
    x = 0;
    y = 0;
}

City::City (double a, double b) {
    x = a;
    y = b;
}

City::~City () {
}

/*
double Old_Distance (City* cityArray, int* order) {
    double d = 0;
    for (int i = 0; i < numOfCities - 1; i++) {
        double dist = (cityArray[order[i+1]].x - cityArray[order[i]].x)*(cityArray[order[i+1]].x - cityArray[order[i]].x);
        dist += (cityArray[order[i+1]].y - cityArray[order[i]].y)*(cityArray[order[i+1]].y - cityArray[order[i]].y);
        d += sqrt(dist);
    }
    double dist = (cityArray[order[numOfCities-1]].x - cityArray[order[0]].x)*(cityArray[order[numOfCities-1]].x - cityArray[order[0]].x);
    dist += (cityArray[order[numOfCities-1]].y - cityArray[order[0]].y)*(cityArray[order[numOfCities-1]].y - cityArray[order[0]].y);
    d += sqrt(dist);
    
    return d;
}*/

double Distance (double* matrix, int* order) {
    double d = 0;
    for (int i = 0; i < numOfCities - 1; i++) {
        d += Dist(matrix, order[i], order[i+1]);
    }
    d += Dist(matrix, order[numOfCities-1], order[0]);
    
    return d;
}

double Distance_2 (double* matrix, int* order) {
    double d = 0;
    for (int i = 0; i < numOfCities - 1; i++) {
        d += Dist(matrix, order[i], order[i+1])*Dist(matrix, order[i], order[i+1]);
    }
    d += Dist(matrix, order[numOfCities-1], order[0])*Dist(matrix, order[numOfCities-1], order[0]);
    
    return d;
}

void Initialize_Distance (double* matrix, City* cityArray) {
    int counter = 0;
    for (int i = 1; i < numOfCities; i++) {
        for (int j = 0; j < i; j++) {
            matrix[counter] = Computed_Distance(cityArray[i], cityArray[j]);
            counter++;
        }
    }
}/*
#include <iostream>
using namespace std;*/
double Dist(double* matrix, int i, int j) {
    if (i == j) return 0;
    if (i < j) return Dist(matrix, j, i);
    int counter = 0;
    for (int k = 1; k < i; k++) {
        counter += k;
    }
//    cout << i << "\t" << j << "\t" << counter+j << endl;
    return matrix[counter + j];
}

double Computed_Distance(City A, City B) {
    double dist = (A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y);
//    return dist;
    return sqrt(dist);
}
