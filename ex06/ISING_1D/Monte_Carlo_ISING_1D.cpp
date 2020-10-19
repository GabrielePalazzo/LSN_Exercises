/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main () {
    Input(); //Inizialization
    if (equilibrate) {
        for (int istep = 1; istep <= 1000; ++istep) {
            Move(metro);
            Measure();
        }
    }
    for (int iblk = 1; iblk <= nblk; ++iblk) { //Simulation
        Reset(iblk);   //Reset block averages
        for (int istep = 1; istep <= nstep; ++istep) {
            Move(metro);
            Measure();
            Accumulate(); //Update block averages
        }
        //cout << "CIAO" << endl;
        Averages(iblk);   //Print results for current block
    }
    ConfFinal(); //Write final configuration

    delete[] blockOutputFileName;
    
    return 0;
}


void Input (void) {
    ifstream ReadInput;
    
    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
  
//Read input informations
    ReadInput.open("input.dat");
    
    ReadInput >> temp;
    beta = 1.0/temp;
    cout << "Temperature = " << temp << endl;
    
    ReadInput >> nspin;
    cout << "Number of spins = " << nspin << endl;
    
    ReadInput >> J;
    cout << "Exchange interaction = " << J << endl;

    ReadInput >> h;
    cout << "External field = " << h << endl << endl;
    
    ReadInput >> metro; // if=1 Metropolis else Gibbs

    ReadInput >> nblk;

    ReadInput >> nstep;
    
    ReadInput >> restart;
    
    ReadInput >> equilibrate;
    
    if(metro==1) cout << "The program perform Metropolis moves" << endl;
    else cout << "The program perform Gibbs moves" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();

    

//Prepare arrays for measurements
    iu = 0; //Energy
    ic = 1; //Heat capacity
    im = 2; //Magnetization
    ix = 3; //Magnetic susceptibility
    
    n_props = 4; //Number of observables
    
    blockOutputFileName = new const char*[n_props];
    
    blockOutputFileName[iu] = "output.ene.0";
    blockOutputFileName[ic] = "output.heat.0";
    blockOutputFileName[im] = "output.mag.0";
    blockOutputFileName[ix] = "output.chi.0";
    
    //for (int i = 0; i < 100000; i++) rnd.Rannyu();    // Aggiunta per cambiare la configurazione iniziale !!!!!!
//initial configuration
    if (restart) {
        ifstream read_restart;
        read_restart.open("config.0");
        for (int i = 0; i < nspin; ++i) {
            read_restart >> s[i];
        }
        read_restart.close();
    }
    else {
        for (int i = 0; i < nspin; ++i) {
            if (rnd.Rannyu() >= 0.5) s[i] = 1;
            else s[i] = -1;
        }
    }
    for (int i = 0; i < nspin; ++i) {
        cout << s[i] << endl;
    }
    //for (int i = 0; i < 200000; i++) rnd.Rannyu();    // Aggiunta per cambiare i numeri casuali generati in seguito !!!!!!
  
//Evaluate energy etc. of the initial configuration
    Measure();

//Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

int Flip (int value) {
    // return value * (-1);
    if (value == 1) return -1;
    else return 1;
}

void Move (int metro) {
    int o;
    //double p, energy_old, energy_new, sm;
    //double energy_up, energy_down;
    double gibbs_probability;
    
    double energyOld, energyNew;
    
    for (int i = 0; i < nspin; ++i) {
        //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
        o = (int)(rnd.Rannyu()*nspin);

        if  (metro == 1) { //Metropolis
// INCLUDE YOUR CODE HERE
            energyOld = Boltzmann(s[o], o);
            energyNew = Boltzmann(Flip(s[o]), o);
            if (energyNew < energyOld) {
                s[o] = Flip(s[o]);
                accepted++;
            }
            else
                if (exp(-beta * (energyNew - energyOld)) > rnd.Rannyu()) {
                    s[o] = Flip(s[o]);
                    accepted++;
                }
            attempted++;
        }
        else { //Gibbs sampling
// INCLUDE YOUR CODE HERE
            //energy_old = Boltzmann(-1, o);
            //energy_new = Boltzmann(1, o);
            // calcolo la probabilità che lo spin per la particella o-esima sia +1
            /*if (temp == 0) {    // per T = 0, l'esponente è infinito
                if (energy_new - energy_old) gibbs_probability = 1;
                else gibbs_probability = 0;
            }*/
            //else {
                //gibbs_probability = 1 / (1 + exp(-2 * beta * J * (Pbc(s[o-1]) + Pbc(s[o+1]))));
            //}
            gibbs_probability = exp(-beta * J * (Boltzmann(1, o))) / (exp(-beta * J * (Boltzmann(1, o))) + exp(-beta * J * Boltzmann(-1, o)));
            //gibbs_probability = 1.0 / (1.0 + exp(-2 * beta * (J * (s[Pbc(o-1)] + s[Pbc(o+1)])) + h));
            if (gibbs_probability > rnd.Rannyu()) {
                s[o] = 1;
                accepted++;
            }
            else {
                s[o] = -1;
                accepted++;
            }
            attempted++;
        }
    }
}

double Boltzmann (int sm, int ip) {
    double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
    return ene;
}

void Measure () {
    //int bin;
    double u = 0.0, m = 0.0;
    //double s_tot = 0;
    //double u_2 = 0;

//cycle over spins
    for (int i = 0; i < nspin; ++i) {
        u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
        m += s[i];
        //u_2 += (-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]))*(-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]));
//        cout << s[i] << endl;
// INCLUDE YOUR CODE HERE
    }
    //cout << u*u << endl;
    //u_2 = u*u;
    walker[iu] = u;
    walker[ic] = u * u;
    walker[im] = m;
    walker[ix] = beta * m * m;
// INCLUDE YOUR CODE HERE
}


void Reset(int iblk) { //Reset block averages
   
    if(iblk == 1) {
        for (int i = 0; i < n_props; ++i) {
            glob_av[i] = 0;
            glob_av2[i] = 0;
        }
    }

    for (int i = 0; i < n_props; ++i) {
        blk_av[i] = 0;
        blk_av2[i] = 0;
    }
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}


void Accumulate(void) { //Update block averages

    for (int i = 0; i < n_props; ++i) {
        blk_av[i] = blk_av[i] + walker[i];
        //blk_av2[i] = blk_av2[i] + (walker[i] * walker[i]);
    }
    blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
    //ofstream Ene, Heat, Mag, Chi;
    const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    
    for (int i = 0; i < n_props; i++) {
        
        ofstream dataBlockingWriter;
        dataBlockingWriter.open(blockOutputFileName[i], ios::app);
                
        double observable = 0;
        if (i != 1) observable = blk_av[i]/blk_norm/(double)nspin;
        else observable = beta * beta * ((blk_av[ic]/blk_norm) - blk_av[iu]/blk_norm * blk_av[iu]/blk_norm) / (double)nspin;
        glob_av[i]  += observable;
        glob_av2[i] += observable * observable;
        double uncertainty = Error(glob_av[i]/(double)iblk,glob_av2[i]/(double)iblk,iblk);
        dataBlockingWriter << setw(wd) << iblk <<  setw(wd) << observable << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << uncertainty << endl;
        
        dataBlockingWriter.close();
        
    }
    /*
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    // Qui
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();
    
    cout << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    //cout << setw(wd) << iblk <<  setw(wd) << glob_av2[iu]/(double)iblk - glob_av[iu]/(double)iblk * glob_av[iu]/(double)iblk << endl;
    
    //stima_c = beta * beta * ((blk_av[ic]/blk_norm/(double)nspin) - stima_u * stima_u); //Heat capacity
    //cout << blk_norm << endl;
    //cout << setw(wd) << beta * beta *(blk_av[ic]/blk_norm - blk_av[iu]/blk_norm * blk_av[iu]/blk_norm) / (double)nspin << endl;
    
    Heat.open("output.heat.0",ios::app);
    stima_c = beta * beta *(blk_av[ic]/blk_norm - blk_av[iu]/blk_norm * blk_av[iu]/blk_norm) / (double)nspin; //Heat capacity
    // Qui
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    //cout << setw(wd) << iblk <<  setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << (glob_av[iu]/(double)iblk) <<  setw(wd) << (glob_av[iu]/(double)iblk) * (glob_av[iu]/(double)iblk) << endl;
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();
    
    Mag.open("output.mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();
    
    Chi.open("output.chi.0",ios::app);
    stima_x = blk_av[ix]/blk_norm/(double)nspin; //Magnetic susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();

// INCLUDE YOUR CODE HERE
*/
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
