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
#include <fstream>
#include <string>
#include "random.h"

using namespace std;
 
#define generatedNumbers 1000000

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

    ofstream write;
    write.open("linear.txt");
    for(int i=0; i<generatedNumbers; i++){
        write << setprecision(10) << rnd.Rannyu() << endl;
    }
    write.close();
    
    write.open("importance.txt");
    for(int i=0; i<generatedNumbers; i++){
        write << setprecision(10) << rnd.Importance_cos() << endl;
    }
    write.close();

    write.open("importance_taylor.txt");
    for(int i=0; i<generatedNumbers; i++){
        write << setprecision(10) << rnd.Importance_taylor() << endl;
    }
    write.close();

    write.open("cubic_lattice.txt");
    for(int i=0; i<generatedNumbers; i++){
        write << setprecision(10) << rnd.Rannyu() << endl;
    }
    write.close();
    
    write.open("continuum.txt");
    for(int i=0; i<2*generatedNumbers; i++){
        write << setprecision(10) << rnd.Rannyu() << endl;
    }
    write.close();

   rnd.SaveSeed();
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
