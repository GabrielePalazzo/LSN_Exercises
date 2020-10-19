#include <iostream>
#include <cmath>
#include <fstream>
#include <cfloat>
#include <ctime>

#include "random.h"
#include "city.h"
#include "genetic.h"

using namespace std;


int seed[4];
Random rnd;

double bestDistance = DBL_MAX;

void Shuffle(int*, int);
void Quick_Sort(double*, int, int);
void Swap(double*, int, int);

int main (int argc, char** argv) {
    
    clock_t t;
    t = clock();
    if (numOfCities < 2) return 1;
    
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
    // (1)
    
    City cities[numOfCities];
    int order[numOfCities];
    int bestOrder[numOfCities];
    int** population = new int*[numOfPopulations];
    for (int i = 0; i < numOfPopulations; i++) {
        population[i] = new int[numOfCities];
    }
    double* distanceMatrix = new double[numOfCities*(numOfCities-1)/2];
    double fitness[numOfPopulations];
    //srand(t);
    
    /*
    int number = rand();
    for (int i = 0; i < number; i++) {
        cities[0].x = rnd.Rannyu();
    }*/
    
    for (int i = 0; i < numOfCities; i++) {
        
#ifdef _BOX_
        cities[i].x = (rnd.Rannyu() - 0.5) * 2.0;
        cities[i].y = (rnd.Rannyu() - 0.5) * 2.0;
#endif
#ifdef _CIRCLE_
        double angle = rnd.Rannyu()*2*M_PI;
        cities[i].x = cos(angle);
        cities[i].y = sin(angle);
#endif
        
        order[i] = i;
    }
    
#ifdef _FILE_
    ifstream read;
    read.open("cities.txt");
    for (int i = 0; i < numOfCities; i++) {
        double num;
        char c[2];
        
        read >> c;
        read >> num;
        cities[i].x = num;
        read >> num;
        cities[i].y = num;
        //cout << cities[i].x << "\t" << cities[i].y << endl;
        order[i] = i;
    }
    read.close();
#endif
    
    //(2)
    
    Initialize_Distance(distanceMatrix, cities);
    
    // (4)
    
    for (int i = 0; i < numOfPopulations; i++) {
        for (int j = 0; j < numOfCities; j++) {
            (population[i])[j] = order[j];
        }
        Shuffle(population[i], 20000);
//        fitness[i] = 1 / (Distance(cities, population[i]) + 1);
//        cout << endl;
        for (int j = 0; j < numOfCities; j++) {
//            cout << order[j] << "\t" << (population[i])[j] << endl;
        }
    }
    
    for (int i = 0; i < numOfPopulations; i++) {
        if (Bad_Individual(population[i])) cout << "Bad individual" << endl;
    }
    
    //(5)
    
    Calculate_Fitness(population, fitness, distanceMatrix);
    
    Normalize_Fitness(fitness);
    
    ofstream writeBestConfiguration, writeAverages;
    
    // (6)
    
    for (int i = 0; i < numOfPopulations; i++) {
        if (Distance(distanceMatrix, population[i]) < bestDistance) {
            bestDistance = Distance(distanceMatrix, population[i]);
            for (int j = 0; j < numOfCities; j++) {
                bestOrder[j] = (population[i])[j];
            }
        }
    }
    
    
    writeAverages.open("path_average.txt");
    
    double* distanceArray = new double[numOfPopulations];
    for (int i = 0; i < numOfPopulations; i++) {
        distanceArray[i] = Distance(distanceMatrix, population[i]);
    }/*
    cout << "----------------------------------------------------------" << endl;
    for (int i = 0; i < numOfPopulations; i++) {
        cout << distanceArray[i] << endl;
    }*/
    cout << "----------------------------------------------------------" << endl;
    
    // (7)
    
    Quick_Sort(distanceArray, 0, numOfPopulations - 1);
    
    double somma = 0;
    double somma_2 = 0;
    for (int i = numOfPopulations / 2; i < numOfPopulations; i++) {
        somma += distanceArray[i];
        somma_2 += distanceArray[i]*distanceArray[i];
    }
    somma /= numOfPopulations - numOfPopulations / 2;
    somma_2 /= numOfPopulations - numOfPopulations / 2;
    
    writeAverages << somma << " " << bestDistance << endl;
    //writeAverages << somma << " " << sqrt(somma_2/numOfPopulations - somma/numOfPopulations*somma/numOfPopulations) << endl;
    
    
    // (8)
    
    Next_Generation(cities, population, fitness, distanceMatrix, &rnd);
    
    for (int i = 0; i < numOfPopulations; i++) {
        if (Bad_Individual(population[i])) cout << "Bad individual" << endl;
    }
    
    for (int j = 0; j < numOfGenerations; j++) {
        Calculate_Fitness(population, fitness, distanceMatrix);
        Normalize_Fitness(fitness);
        Next_Generation(cities, population, fitness, distanceMatrix, &rnd);
        for (int i = 0; i < numOfPopulations; i++) {
            if (Bad_Individual(population[i])) cout << "Bad individual" << endl;
        }
        bool isTrue = false;
        double sum = 0, sum2 = 0;
        for (int i = 0; i < numOfPopulations; i++) {
//            if (Old_Distance(cities, population[i]) < bestDistance) bestDistance = Old_Distance(cities, population[i]);
            sum += Distance(distanceMatrix, population[i]);
            sum2 += Distance(distanceMatrix, population[i])*Distance(distanceMatrix, population[i]);
            if (Distance(distanceMatrix, population[i]) < bestDistance) {
                bestDistance = Distance(distanceMatrix, population[i]);
                isTrue = true;
                for (int j = 0; j < numOfCities; j++) {
                    bestOrder[j] = (population[i])[j];
                }
                cout << j << "\t" << bestDistance << endl;
            }
        }
        if (isTrue || j == numOfGenerations-1) {cout << sum/numOfPopulations << " ± "
            << sqrt(sum2/numOfPopulations - sum/numOfPopulations*sum/numOfPopulations) << endl;
            
//            cout << "\t" << bestDistance << endl;
        }
        
        double currentBestDistance = 100000;
        for (int i = 0; i < numOfPopulations; i++) {
            distanceArray[i] = Distance(distanceMatrix, population[i]);
            if (distanceArray[i] < currentBestDistance) currentBestDistance = distanceArray[i];
        }
        somma = 0;
        somma_2 = 0;
        for (int i = numOfPopulations / 2; i < numOfPopulations; i++) {
            somma += distanceArray[i];
            somma_2 += distanceArray[i]*distanceArray[i];
        }
        somma /= numOfPopulations - numOfPopulations / 2;
        somma_2 /= numOfPopulations - numOfPopulations / 2;
        writeAverages << somma << " " << bestDistance << endl;
        //writeAverages << somma << " " << sqrt(somma_2/numOfPopulations - somma/numOfPopulations*somma/numOfPopulations) << endl;
        
    }
    
    writeAverages.close();
    
    cout << "{" << bestOrder[0];
    for (int k = 1; k < numOfCities; k++) {
        cout <<  ", " << bestOrder[k];
    }
    cout << "};" << endl;
    if (absoluteBest > bestDistance) cout << "Yeeeeee" << endl;
    else  cout << "Noooooo" << endl;
    
    writeBestConfiguration.open("best_configuration.txt");
    for (int i = 0; i < numOfCities; i++) {
        writeBestConfiguration << cities[bestOrder[i]].x << " " << cities[bestOrder[i]].y << endl;
    }
    writeBestConfiguration << cities[bestOrder[0]].x << " " << cities[bestOrder[0]].y << endl;
    writeBestConfiguration.close();
    
    Calculate_Fitness(population, fitness, distanceMatrix);
    cout << endl;
    Normalize_Fitness(fitness);
    
    for (int i = 0; i < numOfPopulations; i++) {
        if (Distance(distanceMatrix, population[i]) < bestDistance) bestDistance = Distance(distanceMatrix, population[i]);
    }
    
    ofstream writeAbsoluteBest;
    writeAbsoluteBest.open("best_distance.txt");
    writeAbsoluteBest << bestDistance << endl;
    writeAbsoluteBest.close();
    
    
    t = clock() - t;
    cout << "Execution time: " << (double)t / CLOCKS_PER_SEC << " s" << endl;
    
    
    delete[] distanceArray;
    
    for (int i = 0; i < numOfPopulations; i++) {
        delete[] population[i];
    }
    delete[] population;
    delete[] distanceMatrix;
    
    return 0;
}

void Shuffle (int* orderArray, int numOfTimes) {
    for (int i = 0; i < numOfTimes; i++) {
        int indexA = (int) rnd.Rannyu(1, numOfCities);
        int indexB = (int) rnd.Rannyu(1, numOfCities);
        Swap(orderArray, indexA, indexB);
    }
}

void Quick_Sort (double* distanceArray, int a, int b) {
    if ((b - a) == 0) return;
    else if (b < a) return Quick_Sort(distanceArray, b, a);
    else {
        int pivotIndex = b;// Scelgo un elemento pivot;
        for (int i = a; i <= b; i++) {
            if ((distanceArray[i] > distanceArray[pivotIndex]) && (pivotIndex > i)) {
                //cout << i << "\t" << distanceArray[i] << "\t" << pivotIndex << "\t" << distanceArray[pivotIndex] << endl;
                Swap(distanceArray, i, pivotIndex);
                pivotIndex = i;
                //cout << i << "\t" << distanceArray[i] << "\t" << pivotIndex << "\t" << distanceArray[pivotIndex] << endl;
            }
            else if ((distanceArray[i] < distanceArray[pivotIndex]) && (pivotIndex < i)) {
                //cout << i << "\t" << distanceArray[i] << "\t" << pivotIndex << "\t" << distanceArray[pivotIndex] << endl;
                Swap(distanceArray, i, pivotIndex);
                int j = pivotIndex;
                pivotIndex = i;
                i = j;
                //cout << i << "\t" << distanceArray[i] << "\t" << pivotIndex << "\t" << distanceArray[pivotIndex] << endl;
            }
        }
        //cout << pivotIndex << endl;
        if (pivotIndex > a) Quick_Sort(distanceArray, a, pivotIndex - 1);
        if (pivotIndex < b) Quick_Sort(distanceArray, pivotIndex + 1, b);
    }
}

void Swap (double* array, int i, int j) {
    double C = array[i];
    array[i] = array[j];
    array[j] = C;
}
