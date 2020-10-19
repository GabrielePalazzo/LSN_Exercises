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

#define generatedNumbersDirect 100000
#define generatedNumbersDiscretized 100000000

//#define generatedNumbersDirect 1000000
//#define generatedNumbersDiscretized 100000000
 
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
    write.open("gauss_t_call.txt");
    for(int i = 0; i < generatedNumbersDirect; i++){
        write << setprecision(10) << rnd.Gauss(0,1) << endl;
    }
    write.close();
    
    write.open("gauss_1_call.txt");
    for(int i = 0; i < generatedNumbersDiscretized; i++){
        write << setprecision(10) << rnd.Gauss(0,1) << endl;
    }
    write.close();
    
    write.open("gauss_t_put.txt");
    for(int i = 0; i < generatedNumbersDirect; i++){
        write << setprecision(10) << rnd.Gauss(0,1) << endl;
    }
    write.close();
    
    write.open("gauss_1_put.txt");
    for(int i = 0; i < generatedNumbersDiscretized; i++){
        write << setprecision(10) << rnd.Gauss(0,1) << endl;
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
