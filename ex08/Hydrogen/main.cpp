/****************************************************************
*****************************************************************
    _/     _/  _/_/_/ _/       Numerical Simulation Laboratory
   _/_/  _/ _/        _/       Physics Department
  _/  _/_/    _/     _/       Universita' degli Studi di Milano
 _/     _/       _/  _/       Prof. D.E. Galli
_/     _/  _/_/_/ _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <iomanip>
//#include <fstream>
#include <string>
#include "random.h"
#include "metropolis.h"
#include "vmc.h"

using namespace std;

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    int accepted = 0;
    double x;
    WaveFunction* psi;
    psi = new Variational(MU, SIGMA);
    ofstream writeConfiguration, writeEnergy;
    writeConfiguration.open("final_configuration.txt");
    writeEnergy.open("energy.txt");
    
    if (argc == 1) {
        cout << "No starting point set. Setting starting point to (0)" << endl;
        x = 0;
    }
    else if (argc == 2) {
        cout << "Setting starting point to (" << atof(argv[1]) << ")" << endl;
        x = atof(argv[1]);
    }
    else {
        cout << "Invalid input. Setting starting point to (0)" << endl;
        x = 0;
    }
    
    //Write(x, writeConfiguration);
    
    //for (int i = 0; i < numberDimensions; i++) cout << x[i] << endl;
    
    for (int i = 0; i < equilibriumSteps; i++) {
        Metropolis_Step(x, rnd, samplingType, psi, accepted);
    }
    accepted = 0;
    //int numberOfBlocks = 200;
    if (MEASURE) {
        int stepsPerBlock = numberOfSteps / numberOfBlocks;
        double averageEnergy = 0;
        double averageEnergySquared = 0;
        for (int i = 0; i < numberOfBlocks; i++) {
            double energyInsideBlock = 0;
            for (int j = 0; j < stepsPerBlock; j++) {
                Metropolis_Step(x, rnd, samplingType, psi, accepted);
                Write(x, writeConfiguration);
                energyInsideBlock += psi->Energy(x);
            }
            energyInsideBlock /= stepsPerBlock;
            averageEnergy += energyInsideBlock;
            averageEnergySquared += energyInsideBlock * energyInsideBlock;
            Write(averageEnergy/(i+1),
                  sqrt((averageEnergySquared/(i+1) - averageEnergy/(i+1) * averageEnergy/(i+1)) / i), writeEnergy);
        }
    }
    if (MINIMIZE) {
        ofstream writeParameters;
        double* grad = new double[2];
        double* param = new double[2];
        grad[0] = 0;
        grad[1] = 0;
        param[0] = MU;
        param[1] = SIGMA;
        double energy = 1000.0;
        for (int i = 0; i < 8000; i++) {
            Minimize(psi, x, rnd, grad, i, param, energy);
        }
        //cout << param[0] << " " << param[1] << " " << energy << endl;
        //cout << param[0] << " " << param[1] << endl;
        writeParameters.open("minimal_parameters.txt");
        //writeParameters << param[0] << " " << param[1] << endl;
        writeParameters << psi->Get_mu() << " " << psi->Get_sigma() << endl;
        writeParameters.close();
        delete[] grad;
        delete[] param;
        //cout << psi->Get_mu() << " " << psi->Get_sigma() << endl;
    }
    cout << "Acceptance rate: " << accepted/(double)numberOfSteps << endl;
    writeConfiguration.close();
    writeEnergy.close();
    delete psi;
     
/*
   rnd.SaveSeed();*/
   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
