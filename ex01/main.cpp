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
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

#define generatedNumbers 10000000

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
    
    write.open("exponential.txt");
    for(int i=0; i<generatedNumbers; i++){
        write << setprecision(10) << rnd.Exponential(1) << endl;
    }
    write.close();
    
    write.open("lorentzian.txt");
    for(int i=0; i<generatedNumbers; i++){
        write << setprecision(10) << rnd.Lorentzian(0,1) << endl;
    }
    write.close();
    
    write.open("buffon_theta.txt");
    for(int i=0; i<generatedNumbers; i++){
        write << setprecision(10) << rnd.Rannyu() << endl;
    }
    write.close();
    
    write.open("buffon_mc.txt");
    for(int i=0; i<generatedNumbers/2; i++){
        write << setprecision(10) << rnd.Rannyu() << endl;
        double x = rnd.Rannyu();
        double y = rnd.Rannyu();
        while ((x*x + y*y) >= 1) {
            x = rnd.Rannyu();
            y = rnd.Rannyu();
        }
        write << setprecision(10) << atan(y/x) << endl;
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
