#include "genetic.h"
#include "random.h"
#include <iostream>
using namespace std;

void Swap (int* orderArray, int indexA, int indexB) {
    //if (indexA == 0 || indexB == 0) return;
    int tempOrder = orderArray[indexA];
    orderArray[indexA] = orderArray[indexB];
    orderArray[indexB] = tempOrder;
}

void Shift (int* orderArray, int numberOfPositionsShifted) {
    
    int* tempCities = new int[numberOfPositionsShifted];
    for (int i = 0; i < numberOfPositionsShifted; i++) {
        tempCities[i] = orderArray[i + 1];
    }
    
    for (int i = 1; i < numOfCities - numberOfPositionsShifted; i++) {
        orderArray[i] = orderArray[i + numberOfPositionsShifted];
    }
    
    for (int i = 0; i < numberOfPositionsShifted; i++) {
        orderArray[numOfCities - numberOfPositionsShifted + i] = tempCities[i];
    }
    
    delete[] tempCities;
}

void Permutation (int* orderArray, int numberOfPositionsShifted, int positionOfFirstGroup, int positionOfSecondGroup) {
    int* tempCities = new int[numberOfPositionsShifted];
    for (int i = 0; i < numberOfPositionsShifted; i++) {
        tempCities[i] = orderArray[i + positionOfFirstGroup];
    }
    
    for (int i = 0; i < numberOfPositionsShifted; i++) {
        orderArray[i + positionOfFirstGroup] = orderArray[i + positionOfSecondGroup];
    }
    
    for (int i = 0; i < numberOfPositionsShifted; i++) {
        orderArray[i+positionOfSecondGroup] = tempCities[i];
    }
    delete[] tempCities;
}

void Inversion (int* orderArray, int numberOfPositionsShifted, int positionOfGroup) {
    int* tempCities = new int[numberOfPositionsShifted];
    for (int i = 0; i < numberOfPositionsShifted; i++) {
        tempCities[i] = orderArray[i + positionOfGroup];
    }
    
    for (int i = 0; i < numberOfPositionsShifted; i++) {
        orderArray[i + positionOfGroup] = tempCities[numberOfPositionsShifted - i - 1];
    }
    
    delete[] tempCities;
}

void Smart_Swap (City* cityArray, int* pickedOrder, Random* rnd) {
    bool found = false;
    double xa, xb, xc, xd, ya, yb, yc, yd;
    int indexA = (int) rnd->Rannyu(1, numOfCities);
    int indexB = (int) rnd->Rannyu(1, numOfCities);
    //int indexA = (int) Rand_uniform(1, numOfCities);
    //int indexB = (int) Rand_uniform(1, numOfCities);
    for (int i = 0; i < numOfCities - 2 && !found; i++) {
        xa = cityArray[pickedOrder[i]].x;
        ya = cityArray[pickedOrder[i]].y;
        xb = cityArray[pickedOrder[i+1]].x;
        yb = cityArray[pickedOrder[i+1]].y;
        for (int j = i + 1; j < numOfCities - 1 && !found; j++) {
            xc = cityArray[pickedOrder[j]].x;
            yc = cityArray[pickedOrder[j]].y;
            xd = cityArray[pickedOrder[j+1]].x;
            yd = cityArray[pickedOrder[j+1]].y;
            if (Intersect(cityArray, pickedOrder[i], pickedOrder[i+1], pickedOrder[j], pickedOrder[j+1])) {
                found = true;
                indexA = i + 1;
                indexB = j + 1;
                break;
            }
//            cout << Intersect(cityArray, pickedOrder[i], pickedOrder[i+1], pickedOrder[j], pickedOrder[j+1]) << endl;
        }
    }
    Swap(pickedOrder, indexA, indexB);
}

void Random_Swap (int* pickedOrder, Random* rnd) {
    int indexA = (int) rnd->Rannyu(1, numOfCities);
    int indexB = (int) rnd->Rannyu(1, numOfCities);
    //int indexA = (int) Rand_uniform(1, numOfCities);
    //int indexB = (int) Rand_uniform(1, numOfCities);
    while (indexB == indexA) indexB = (int) Rand_uniform(1, numOfCities);
    Swap(pickedOrder, indexA, indexB);
}

void Random_Shift (int* pickedOrder, Random* rnd) {
    int numberOfPositionsShifted = (int) rnd->Rannyu(1, numOfCities - 2);
    //int numberOfPositionsShifted = (int) Rand_uniform(1, numOfCities - 2); // Is N-2 better than N-1?
    Shift(pickedOrder, numberOfPositionsShifted);
}

void Random_Permutation (int* pickedOrder, Random* rnd) {
    int numberOfPositionsShifted = (int) rnd->Rannyu(1, (numOfCities - 1) / 2);
    int positionOfFirstGroup = (int) rnd->Rannyu(1, numOfCities - 2 * numberOfPositionsShifted);
    int positionOfSecondGroup = (int) rnd->Rannyu(positionOfFirstGroup + numberOfPositionsShifted, numOfCities - numberOfPositionsShifted);
    /*int numberOfPositionsShifted = (int) Rand_uniform(1, (numOfCities - 1) / 2);
    int positionOfFirstGroup = (int) Rand_uniform(1, numOfCities - 2 * numberOfPositionsShifted);
    int positionOfSecondGroup = (int) Rand_uniform(positionOfFirstGroup + numberOfPositionsShifted, numOfCities - numberOfPositionsShifted);*/
    Permutation(pickedOrder, numberOfPositionsShifted, positionOfFirstGroup, positionOfSecondGroup);
}

void Random_Inversion (int* pickedOrder, Random* rnd) {
    int numberOfPositionsShifted = (int) rnd->Rannyu(1, numOfCities - 1);
    int positionOfGroup = (int) rnd->Rannyu(1, numOfCities - numberOfPositionsShifted);
    //int numberOfPositionsShifted = (int) Rand_uniform(1, numOfCities - 1);
    //int positionOfGroup = (int) Rand_uniform(1, numOfCities - numberOfPositionsShifted);
    Inversion(pickedOrder, numberOfPositionsShifted, positionOfGroup);
}

void Calculate_Fitness (int** orderArray, double* fit, double* matrix) {
    for (int i = 0; i < numOfPopulations; i++) {
//        fit[i] = 1. / (Old_Distance(cityArray, orderArray[i]) + 1.);
        fit[i] = 1. / (Distance(matrix, orderArray[i]) + 1.);
    }
}

void Normalize_Fitness (double* fit) {
    double sum = 0;
    for (int i = 0; i < numOfPopulations; i++) {
        sum += fit[i];
    }
    
    for (int i = 0; i < numOfPopulations; i++) {
        fit[i] /= sum;
    }
}

void Next_Generation (City* cityArray, int** previousPopulation, double* fit, double* matrix, Random* rnd) {
    int** newPopulation = new int*[numOfPopulations];
    for (int i = 0; i < numOfPopulations; i++) {
        newPopulation[i] = new int[numOfCities];
    }
    
    for (int i = 0; i < numOfPopulations; i++) {
        int pickedIndex = Pick(fit, rnd);
        
        for (int j = 0; j < numOfCities; j++) {
            (newPopulation[i])[j] = (previousPopulation[pickedIndex])[j];
        }
    }
    
    for (int i = 0; i < numOfPopulations; i++) {
        
        //double probability = Rand_uniform();
        double probability = rnd->Rannyu();
        if (probability < crossoverRate && i < numOfPopulations - 1) {
//            int pickedIndexA = Pick(fit);
            int pickedIndexA = Pick_Max(fit, rnd);
            int pickedIndexB = Pick(fit, rnd);
            while (pickedIndexB == pickedIndexA) pickedIndexB = Pick(fit, rnd);
            int* testPopulationA = new int[numOfCities];
            int* testPopulationB = new int[numOfCities];
            int* testPopulation = new int[numOfCities];
            for (int j = 0; j < numOfCities; j++) {
                testPopulationA[j] = (previousPopulation[pickedIndexA])[j];
                testPopulationB[j] = (previousPopulation[pickedIndexB])[j];
            }
            Crossover(newPopulation[i], newPopulation[i+1], testPopulationA, testPopulationB, rnd);
            i++;
            delete[] testPopulationA;
            delete[] testPopulationB;
            delete[] testPopulation;
        }
    }
    
    
//    Pick(previousPopulation, fit);
    double avg_fit = 0, max_fit = 0;
    for (int i = 0; i < numOfPopulations; i++) {
        double new_fit = 1. / (Distance_2(matrix, newPopulation[i]) + 1.);
        avg_fit += new_fit;
        if (new_fit > max_fit) max_fit = new_fit;
    }
    avg_fit /= numOfPopulations;
    if (CHINESE) {
        for (int i = 0; i < numOfPopulations; i++) {
            //double probability = Rand_uniform();
            double probability = rnd->Rannyu();
            double new_fit = 1. / (Distance_2(matrix, newPopulation[i]) + 1.);
            new_fit /= numOfPopulations;
            if (new_fit < avg_fit) {
                if (probability < maxMutationRate) Mutate(newPopulation[i], rnd);
            }
            else {
                if (probability < 1*(maxMutationRate - minMutationRate)*(new_fit-avg_fit)/(max_fit-avg_fit)) Mutate(newPopulation[i], rnd);
            }
        }
    }
    else {
        for (int i = 0; i < numOfPopulations; i++) {
            //double probability = Rand_uniform();
            double probability = rnd->Rannyu();
            if (probability < mutationRate) Mutate(newPopulation[i], rnd);
        }
    }
    
    for (int i = 0; i < numOfPopulations; i++) {
        for (int j = 0; j < numOfCities; j++) {
            (previousPopulation[i])[j] = (newPopulation[i])[j];
        }
    }
    for (int i = 0; i < numOfPopulations; i++) {
        delete[] newPopulation[i];
    }
    delete[] newPopulation;
}

void Simulated_Annealing (int** previousPopulation, double* fit, double* matrix, double T, Random* rnd) {
    int** newPopulation = new int*[numOfPopulations];
    for (int i = 0; i < numOfPopulations; i++) {
        newPopulation[i] = new int[numOfCities];
    }
    for (int i = 0; i < numOfPopulations; i++) {
        for (int j = 0; j < numOfCities; j++) {
            (newPopulation[i])[j] = (previousPopulation[i])[j];
        }
    }
    
    for (int i = 0; i < numOfPopulations; i++) {
        Mutate(newPopulation[i], rnd);
    }
    
    for (int i = 0; i < numOfPopulations; i++) {
        double delta = Distance(matrix, newPopulation[i]) - Distance(matrix, previousPopulation[i]);
        if (delta < 0) {
            for (int j = 0; j < numOfCities; j++) {
                (previousPopulation[i])[j] = (newPopulation[i])[j];
            }
        }
        else {
            //double randomNumber = Rand_uniform();
            double randomNumber = rnd->Rannyu();
            double probability = exp(-delta/T);
            if (randomNumber < probability) {
                for (int j = 0; j < numOfCities; j++) {
                    (previousPopulation[i])[j] = (newPopulation[i])[j];
                }
            }
        }
    }
    for (int i = 0; i < numOfPopulations; i++) {
        delete[] newPopulation[i];
    }
    delete[] newPopulation;
}

int Pick (double* fit, Random* rnd) {
    int index = 0;
    //double probability = Rand_uniform();
    double probability = rnd->Rannyu();
//    cout << probability << endl;
    while (probability > 0) {
        probability -= fit[index];
        index++;
    }
    index--;
    return index;
//    return previousPopulation[index];
}

int Pick_Max (double* fit, Random* rnd) {
    double current_max = 0;
    double previous_max = 0;
    double next_previous_max = 0;
    int index = 0;
    int indexA = 0;
    int indexB = 0;
    int indexC = 0;
    for (int i = 0; i < numOfPopulations; i++) {
        if (fit[i] > current_max) {
            next_previous_max = previous_max;
            indexC = indexB;
            previous_max = current_max;
            indexB = indexA;
            current_max = fit[i];
            indexA = i;
            index = i;
        }
        else if (fit[i] > previous_max) {
            next_previous_max = previous_max;
            indexC = indexB;
            previous_max = fit[i];
            indexB = i;
        }
        else if (fit[i] > next_previous_max) {
            next_previous_max = fit[i];
            indexC = i;
        }
    }return indexA;
    //double probability = Rand_uniform();
    double probability = rnd->Rannyu();
    if (probability > 0.5) return indexA;
    else if (probability > 0.8)  return indexB;
    else return indexC;
    
    return index;
    //    return previousPopulation[index];
}

void Mutate (int* pickedOrder, Random* rnd) {
    
    double probability = rnd->Rannyu();
    if (probability < 0.25) {
        double repetition = rnd->Gauss(0, 1);
        if (repetition < 0) repetition *= -1;
        for (int i = 0; i < repetition + 1; i++)
            Random_Swap(pickedOrder, rnd);
    }
    
    probability = rnd->Rannyu();
    if (probability < 0.5) {
        double repetition = rnd->Gauss(0, 1);
        if (repetition < 0) repetition *= -1;
        for (int i = 0; i < repetition; i++)
            Random_Shift(pickedOrder, rnd);
    }
    
    probability = rnd->Rannyu();
    if (probability < 0.5) {
        double repetition = rnd->Gauss(0, 1);
        if (repetition < 0) repetition *= -1;
        for (int i = 0; i < repetition; i++)
            Random_Permutation(pickedOrder, rnd);
    }
    
    probability = rnd->Rannyu();
    if (probability < 0.5) {
        double repetition = rnd->Gauss(0, 1);
        if (repetition < 0) repetition *= -1;
        for (int i = 0; i < repetition; i++)
            Random_Inversion(pickedOrder, rnd);
    }
}

void Crossover (int* childA, int* childB, int* father, int* mother, Random* rnd) {
//    childA[0] = father[0]; // == mother[0] (and fixed)
//    childB[0] = mother[0]; // == father[0] (and fixed)
    //int cutIndex = (int) Rand_uniform(1, numOfCities-1);
    int cutIndex = (int) rnd->Rannyu(1, numOfCities-1);
    for (int i = 0; i < cutIndex; i++) {
        childA[i] = father[i];
        childB[i] = mother[i];
    }
    
    
    for (int i = cutIndex; i < numOfCities; i++) {
        for (int j = 1; j < numOfCities; j++) {
            bool found = false;
            int k = 1;
            while (k < cutIndex) {
                if (mother[j] == father[k]) {
                    found = true;
                    break;
                }
                k++;
            }
            if (!found) {
                if (i < numOfCities) {
                    childA[i] = mother[j];
                    i++;
                }
            }
        }
    }
    for (int i = cutIndex; i < numOfCities; i++) {
        for (int j = 1; j < numOfCities; j++) {
            bool found = false;
            int k = 1;
            while (k < cutIndex) {
                if (father[j] == mother[k]) {
                    found = true;
                    break;
                }
                k++;
            }
            if (!found) {
                if (i < numOfCities) {
                    childB[i] = father[j];
                    i++;
                }
            }
        }
    }
}

bool Intersect (City* cityArray, int a, int b, int c, int d) {
    
    int dir1 = Direction(cityArray, a, b, c);
    int dir2 = Direction(cityArray, a, b, d);
    int dir3 = Direction(cityArray, c, d, a);
    int dir4 = Direction(cityArray, c, d, b);
    
    if(dir1 != dir2 && dir3 != dir4)
        return true; //they are intersecting
    
    if(dir1==0 && On_Line(cityArray, a, b, c)) //when p2 of line2 are on the line1
        return true;
    
    if(dir2==0 && On_Line(cityArray, a, b, d)) //when p1 of line2 are on the line1
        return true;
    
    if(dir3==0 && On_Line(cityArray, c, d, a)) //when p2 of line1 are on the line2
        return true;
    
    if(dir4==0 && On_Line(cityArray, c, d, b)) //when p1 of line1 are on the line2
        return true;
    
    return false;
}

int Direction (City* cityArray, int a, int b, int c) {
    double value = (cityArray[b].y-cityArray[a].y)*(cityArray[c].x-cityArray[b].x)-(cityArray[b].x-cityArray[a].x)*(cityArray[c].y-cityArray[b].y);
    if (value == 0) return 0;
    else if (value < 0) return -1;
    else return 1;
    // Collinear = 0, Clockwise > 0, Counterclockwise < 0
}

bool On_Line(City* cityArray, int a, int b, int c) {   //check whether p is on the line or not
    if(cityArray[c].x <= max(cityArray[a].x, cityArray[b].x) && cityArray[c].x <= min(cityArray[a].x, cityArray[b].x) &&
       (cityArray[c].y <= max(cityArray[a].y, cityArray[b].y) && cityArray[c].y <= min(cityArray[a].y, cityArray[b].y)))
        return true;
    
    return false;
}

long power(int a, int n){
    long result = 1;
    for (int i=0;i<n;i++) result*=a;
    return result;
}

int Bad_Individual (int* individual) {
    
    int result = 0;
    
    if (individual[0] != 0) return 0;
    
    // Plan A: Sort the array and check if x_i = i (way slower than plan B)
    if (numOfCities > 64) {
        result = numOfCities - 1;
        int* newIndividual = new int[numOfCities - 1];
        for (int i = 1; i < numOfCities; i++) {
            newIndividual[i - 1] = individual[i];
        }
        Sort_Chromosomes(newIndividual, 0, numOfCities - 2);
        for (int i = 1; i < numOfCities; i++) {
            if (newIndividual[i - 1] == i) result--;
        }
        delete[] newIndividual;
    }
    // Plan B: Use some sort of binary encoding x_i => 2^(x_i). Then i compose the number (sum over all chromosomes (cities). If the number must be equal to all ones: 11111...111111
    else {
        long encoded = 0;
        for (int i = 1; i < numOfCities; i++) {
            //encoded += (long)((long)1<<(individual[i]));  // Way faster than power() and pow()
            encoded |= (long)((long)1<<(individual[i]));  // Slightly faster than with += and less risk of overflow
            //encoded += power(2, individual[i]);
            //encoded += pow(2, individual[i]);   // Slightly slower than power()
        }
        result = (encoded - (long)((long)1<<(numOfCities)) + (long)2);
    }
    return result;
}

void Sort_Chromosomes (int* individual, int a, int b) {
    if ((b - a) == 0) return;
    else if (b < a) return Sort_Chromosomes(individual, b, a);
    else {
        int pivotIndex = b;// Scelgo un elemento pivot;
        for (int i = a; i <= b; i++) {
            if ((individual[i] > individual[pivotIndex]) && (pivotIndex > i)) {
                //cout << i << "\t" << individual[i] << "\t" << pivotIndex << "\t" << individual[pivotIndex] << endl;
                Swap_Chromosomes(individual, i, pivotIndex);
                pivotIndex = i;
                //cout << i << "\t" << distanceArray[i] << "\t" << pivotIndex << "\t" << distanceArray[pivotIndex] << endl;
            }
            else if ((individual[i] < individual[pivotIndex]) && (pivotIndex < i)) {
                //cout << i << "\t" << distanceArray[i] << "\t" << pivotIndex << "\t" << distanceArray[pivotIndex] << endl;
                Swap_Chromosomes(individual, i, pivotIndex);
                int j = pivotIndex;
                pivotIndex = i;
                i = j;
                //cout << i << "\t" << distanceArray[i] << "\t" << pivotIndex << "\t" << distanceArray[pivotIndex] << endl;
            }/*
            for (int j = 0; j < numOfCities; j++) cout << individual[j] << " ";
            cout << endl;*/
        }
        //cout << pivotIndex << endl;
        if (pivotIndex > a) Sort_Chromosomes(individual, a, pivotIndex - 1);
        if (pivotIndex < b) Sort_Chromosomes(individual, pivotIndex + 1, b);
    }
}

void Swap_Chromosomes (int* individual, int i, int j) {
    int C = individual[i];
    individual[i] = individual[j];
    individual[j] = C;
}
