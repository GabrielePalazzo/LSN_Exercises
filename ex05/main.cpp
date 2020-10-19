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
    double* x = new double[numberDimensions];
    WaveFunction* psi;
    if (state == GROUND) {
        psi = new Y_1_0_0();
        //cout << "Psi_{1,0,0}" << endl;
        cout << "|1,0,0>" << endl;
    }
    else if (state == EXCITED){
        psi = new Y_2_1_0();
        //cout << "Psi_{2,1,0}" << endl;
        cout << "|2,1,0>" << endl;
    }
    ofstream writeConfiguration;
    writeConfiguration.open("final_configuration.txt");
    
    if (argc == 1) {
        cout << "No starting point set. Setting starting point to (0, 0, 0)" << endl;
        for (int i = 0; i < numberDimensions; i++) x[i] = 0;
    }
    else if (argc == 4) {
        cout << "Setting starting point to (";
        for (int i = 0; i < numberDimensions; i++) {
            x[i] = atof(argv[i + 1]);
            if (i < numberDimensions - 1) cout << atof(argv[i + 1]) << ", ";
            else cout << atof(argv[i + 1]) << ")" << endl;
        }
    }
    else {
        cout << "Invalid input. Setting starting point to (0, 0, 0)" << endl;
        for (int i = 0; i < numberDimensions; i++) x[i] = 0;
    }
    
    //Write(x, writeConfiguration);
    
    //for (int i = 0; i < numberDimensions; i++) cout << x[i] << endl;
    
    for (int i = 0; i < equilibriumSteps; i++) {
        Metropolis_Step(x, rnd, samplingType, psi, accepted);
    }
    accepted = 0;
    for (int i = 0; i < numberOfSteps; i++) {
        Metropolis_Step(x, rnd, samplingType, psi, accepted);
        Write(x, writeConfiguration);
    }
    
    cout << "Acceptance rate: " << accepted/(double)numberOfSteps << endl;
    writeConfiguration.close();
    delete psi;
    delete[] x;
    
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
