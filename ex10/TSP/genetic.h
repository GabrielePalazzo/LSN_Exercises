#ifndef __genetic_h__
#define __genetic_h__

#include "city.h"
#include "random.h"

void Swap(int*, int, int);
void Shift(int*, int);
void Permutation(int*, int, int, int);
void Inversion(int*, int, int);
void Random_Swap(int*, Random*);
void Smart_Swap(City*, int*, Random*);
void Random_Shift(int*, Random*);
void Random_Permutation(int*, Random*);
void Random_Inversion(int*, Random*);
void Calculate_Fitness(int**, double*, double*);
void Normalize_Fitness(double*);
void Next_Generation(City*, int**, double*, double*, Random*);
void Simulated_Annealing(int**, double*, double*, double, Random*);
int Pick(double*, Random*);
int Pick_Max(double*, Random*);
void Mutate(int*, Random*);
void Crossover(int*, int*, int*, int*, Random*);
bool Intersect(City*, int, int, int, int);
int Direction(City*, int, int, int);
bool On_Line(City*, int, int, int);
int Bad_Individual(int*);
void Sort_Chromosomes(int*, int, int);
void Swap_Chromosomes(int*, int, int);

#endif

