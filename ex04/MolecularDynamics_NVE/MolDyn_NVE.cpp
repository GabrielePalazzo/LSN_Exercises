/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
    Input();             //Inizialization
    int nconf = 1;
    
    //nstep /= blockNumber;
    for (int iblk = 1; iblk <= blockNumber; ++iblk) { //Simulation
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep){
            
            Move();           //Move particles with Verlet algorithm
            if((istep + iblk * nstep)%iprint == 0) cout << "Number of time-steps: " << (istep + (iblk - 1) * nstep) << endl;
            if ((istep + iblk * nstep)% printStep == 0){
                //cout << "Writing " << istep/10 << endl;
                Measure();     //Properties measurement
                Accumulate(); //Update block averages
                
                //cout << stima_pot << endl;
                
                //            ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                nconf += 1;
            }
            
        }
        Averages(iblk);   //Print results for current block
    }
    ConfFinal();         //Write final configuration to restart
    
    delete[] blockOutputFileName;
    
    return 0;
}


void Input(void){ //Prepare all stuff for the simulation
    ifstream ReadInput,ReadConf;
    //double ep, ek, pr, et, vir;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

    seed = 1;    //Set seed for random numbers
    srand(seed); //Initialize random number generator
  
    ReadInput.open("input.dat"); //Read input

    ReadInput >> temp;

    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol,1.0/3.0);
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;
    ReadInput >> readOld;
    ReadInput >> printStep;
    ReadInput >> blockNumber;
    ReadInput >> dataBlocking;
    
    if (dataBlocking == 0) {
        blockNumber = 1;
        cout << "No blocking" << endl;
    }
    
    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl;
    ReadInput.close();
    cout << "Read old configuration from file : ";
    if (readOld) cout << "Yes" << endl;
    else cout << "No" << endl;
    cout << endl;
    ReadInput.close();

//Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    ip = 4; //Pressure
    n_props = 5; //Number of observables
    
    blockOutputFileName = new const char*[n_props];
    
    blockOutputFileName[iv] = "ave_epot.out";
    blockOutputFileName[ik] = "ave_ekin.out";
    blockOutputFileName[ie] = "ave_etot.out";
    blockOutputFileName[it] = "ave_temp.out";
    blockOutputFileName[ip] = "ave_pres.out";
    
    bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
    // (1)
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
    }
    ReadConf.close();
    
    
    ReadConf.open("old.0");
    if (ReadConf.fail()) {
// If old.0 is missing, generate a "new" old configuration
        readOld = false;
        ReadConf.close();
    
    }

// Read old configuration
    if (readOld) {
        // (2)
        cout << "Read old spatial configuration from file old.0 " << endl << endl;
        for (int i=0; i<npart; ++i){
            
            ReadConf >> xold[i] >> yold[i] >> zold[i];
            
            xold[i] = xold[i] * box;
            yold[i] = yold[i] * box;
            zold[i] = zold[i] * box;
            
        }
        
        double newCoordinate, fx[m_part], fy[m_part], fz[m_part];

        for(int i=0; i<npart; ++i){ //Force acting on particle i
            
            fx[i] = Force(i,0);
            fy[i] = Force(i,1);
            fz[i] = Force(i,2);
            
        }

        double sumv2 = 0.0, fs;
        
        for(int i=0; i<npart; ++i){
            
            // Before:
            //   - xold = r[t-dt]
            //   - x = r[t]
            // After:
            //   - xold = r[t]
            //   - x = r[t+dt]
            
// Compute the new coordinates from r[t] and r[t-dt] with Verlet algorithm
            // (3)
//            newCoordinate = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
            newCoordinate = Verlet_step(x[i], xold[i], fx[i], delta);
            xold[i] = x[i];
            x[i] = newCoordinate;
            newCoordinate = Verlet_step(y[i], yold[i], fy[i], delta);
//            newCoordinate = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
            yold[i] = y[i];
            y[i] = newCoordinate;
            newCoordinate = Verlet_step(z[i], zold[i], fz[i], delta);
//            newCoordinate = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
            zold[i] = z[i];
            z[i] = newCoordinate;
            
// Compute the velocities at time t+(dt/2)
            
            vx[i] = Pbc(x[i] - xold[i]) / (delta);
            vy[i] = Pbc(y[i] - yold[i]) / (delta);
            vz[i] = Pbc(z[i] - zold[i]) / (delta);
            
            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
            
        }
        // (4)
        sumv2 /= (double)npart;
//        cout << temp << "\t" << sumv2/3 << endl;
        // (5)
        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        
//        sumv2 = 0;
        for(int i=0; i<npart; ++i){

// Rescale the velocities to match the desired temperature
            
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;
            
//            sumv2 += velx[i]*velx[i] + vely[i]*vely[i] + velz[i]*velz[i];
            
        }
//        sumv2 /= (double)npart;
//        cout << temp << "\t" << sumv2/3 << endl;
        
        for (int i=0; i<npart; ++i){
            // (6)
// Rescale the old coordinates (at time t)
            
            xold[i] = Pbc(x[i] - delta*vx[i]);
            yold[i] = Pbc(y[i] - delta*vy[i]);
            zold[i] = Pbc(z[i] - delta*vz[i]);
            
        }
        
        
        sumv2 = 0;
        for(int i=0; i<npart; ++i){
            
// Check the temperature
            
            vx[i] = Pbc(x[i] - xold[i]) / (delta);
            vy[i] = Pbc(y[i] - yold[i]) / (delta);
            vz[i] = Pbc(z[i] - zold[i]) / (delta);
            
            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
            
        }
        
        sumv2 /= (double)npart;
        //cout << temp << "\t" << sumv2/3 << endl;
        
        
    }
    else {
//Prepare initial velocities (if i haven't read old.0)
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
   
        for (int i=0; i<npart; ++i){
        
            vx[i] = rand()/double(RAND_MAX) - 0.5;
            vy[i] = rand()/double(RAND_MAX) - 0.5;
            vz[i] = rand()/double(RAND_MAX) - 0.5;

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
    
        }
    
        for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    
        double sumv2 = 0.0, fs;
        
        for (int i=0; i<npart; ++i){
        
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];

            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        
        }
        
        sumv2 /= (double)npart;
        //cout << temp << "\t" << sumv2/3 << "\tError" << endl;
        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        
        for (int i=0; i<npart; ++i){
            
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;

            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
            
        }
        sumv2 =0;
        for(int i=0; i<npart; ++i){
            
            // Check the temperature

            vx[i] = Pbc(x[i] - xold[i]) / (delta);
            vy[i] = Pbc(y[i] - yold[i]) / (delta);
            vz[i] = Pbc(z[i] - zold[i]) / (delta);
  
            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
            
        }
        
        sumv2 /= (double)npart;
        //cout << temp << "\t" << sumv2/3 << endl;
        
    }
    //cout << temp << endl;
    
    ReadConf.close();
    
    return;
}


void Move(void){ //Move particles with Verlet algorithm
    double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i
        fx[i] = Force(i,0);
        fy[i] = Force(i,1);
        fz[i] = Force(i,2);
    }
    double summm = 0;

    for(int i=0; i<npart; ++i){ //Verlet integration scheme

        xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
        ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
        znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

        vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
        vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
        vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

        xold[i] = x[i];
        yold[i] = y[i];
        zold[i] = z[i];

        x[i] = xnew;
        y[i] = ynew;
        z[i] = znew;
        
        summm += (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])/npart;
    }
    //cout << temp << "\t" << summm/3 << endl;
    return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
//void Measure(double energyToSave){ //Properties measurement
    int bin;
    double v, t, vij;
    double dx, dy, dz, dr;
    ofstream Epot, Ekin, Etot, Temp, Pres, Gave;
    
    Epot.open("output_epot.dat",ios::app);
    Ekin.open("output_ekin.dat",ios::app);
    Temp.open("output_temp.dat",ios::app);
    Etot.open("output_etot.dat",ios::app);
    Pres.open("output_pres.dat",ios::app);
    //Gave.open("output_gave.dat",ios::app);
    
    v = 0.0; //reset observables
    t = 0.0;
    stima_pres = 0.0;
    
    double wij, w=0.0;
    
    double half_box = box/2.0;
    double half_box_2 = half_box*half_box;
    double r_2;
    
    double r_max, r_min;
    
    for (int k = 0; k < nbins; ++k) histogram[k] = 0.0;

//cycle over pairs of particles
    for (int i=0; i<npart-1; ++i){
        for (int j=i+1; j<npart; ++j){
            
            dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
            dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
            dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)
            
            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);
            
            if (dx > box/2.0) dx -= half_box;
            else if (dx < -box/2.0) dx += half_box;
            if (dy > box/2.0) dy -= half_box;
            else if (dy < -box/2.0) dy += half_box;
            if (dz > box/2.0) dz -= half_box;
            else if (dz < -box/2.0) dz += half_box;
            r_2 = dx*dx + dy*dy + dz*dz;
            if (r_2 < half_box_2) {
                bin = sqrt(r_2)/bin_size;
                r_min = bin*bin_size;
                r_max = (bin+1)*bin_size;
                double rrr = r_max*r_max*r_max - r_min*r_min*r_min;
                histogram[bin] += 2.0*3/(4*M_PI*rrr);
            }
            
            if(dr < rcut){
                vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
                
                wij = (16.0/pow(dr,12) - 8.0/pow(dr,6));

//Potential energy
                v += vij;
                
                w += wij;
            }
        }
    }

//Kinetic energy
    for (int i=0; i<npart; ++i)
        t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle
    
// Pressure
    
    /*for (int i = 0; i < npart-1; i++) {
        for (int j = i+1; j < npart; j++) {
            
            dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
            dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
            dz = Pbc( zold[i] - zold[j] );
            
            dr = dx*dx + dy*dy + dz*dz;
            //dr = sqrt(dr);
            
            if(dr < rcut*rcut){
                //stima_pres += ( 1. / pow((dr), 3) - 0.5) / pow((dr), 3);
                //stima_pres += (16.0/pow(sqrt(dr),12) - 8.0/pow(sqrt(dr),6));
                
            }
                
            
        }
        stima_pres /= (npart - i - 1);
    }*/
    
    //stima_pres /= (npart);
    stima_pres = w;
    stima_pres *= 1.0 / vol;
    
    //cout << stima_pres << "\t";
    
    stima_pres += rho * stima_temp;
    stima_pres = rho * temp + w / vol;
    
    //cout << stima_pres << endl;
    
    //cout << temp << "\t" << stima_temp << endl;
    //cout << stima_etot << "\t" << stima_kin+stima_pot << endl;
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl;
    
    
    walker[iv] = stima_pot;
    walker[ik] = stima_kin;
    walker[it] = stima_temp;
    walker[ie] = stima_etot;
    walker[ip] = stima_pres;
    /*
    for (int i = 0; i < nbins; i++) {
        stima_gdir = histogram[i]/rho/(double)npart;
        //glob_av[i] += stima_gdir;
        //glob_av2[i] += stima_gdir*stima_gdir;
        //err_gdir = Error(glob_av[i],glob_av2[i],iblk);
        Gave << stima_gdir << endl;
        //Gave << setw(wd) << iblk << setw(wd) << stima_gdir << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << err_gdir << endl;
    }
    */
    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();
    //Gave.close();
    
    return;
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
        //blk_av2[i] = 0;
    }
    blk_norm = 0;
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
    const int wd=12;
    cout << "Block number " << iblk << endl;
    
    for (int i = 0; i < n_props; i++) {
        
        ofstream dataBlockingWriter;
        dataBlockingWriter.open(blockOutputFileName[i], ios::app);
        
        double observable = blk_av[i]/blk_norm;
        glob_av[i]  += observable;
        glob_av2[i] += observable * observable;
        double uncertainty = Statistical_uncertainty(glob_av[i]/(double)iblk,glob_av2[i]/(double)iblk,iblk);
        dataBlockingWriter << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << uncertainty << endl;
        
        dataBlockingWriter.close();
        
    }
    
    /*
    ofstream Epot_w, Ekin_w, Etot_w, Temp_w, Pres_w;
    
    Epot_w.open("ave_epot.out", ios::app);
    Ekin_w.open("ave_ekin.out", ios::app);
    Temp_w.open("ave_temp.out", ios::app);
    Etot_w.open("ave_etot.out", ios::app);
    Pres_w.open("ave_pres.out", ios::app);
    
    const int wd=12;
    cout << "Block number " << iblk << endl;
    
    double stima_v = blk_av[iv]/blk_norm; // Potential Energy
    glob_av[iv]  += stima_v;
    glob_av2[iv] += stima_v * stima_v;
    err_v = Statistical_uncertainty(glob_av[iv]/(double)iblk,glob_av2[iv]/(double)iblk,iblk);
    Epot_w << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_v << endl;
    
    double stima_k = blk_av[ik]/blk_norm; // Kinetic Energy
    glob_av[ik]  += stima_k;
    glob_av2[ik] += stima_k * stima_k;
    err_k = Statistical_uncertainty(glob_av[ik]/(double)iblk,glob_av2[ik]/(double)iblk,iblk);
    Ekin_w << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_k << endl;
    
    double stima_t = blk_av[it]/blk_norm; // Temperature
    glob_av[it]  += stima_t;
    glob_av2[it] += stima_t * stima_t;
    err_t = Statistical_uncertainty(glob_av[it]/(double)iblk,glob_av2[it]/(double)iblk,iblk);
    Temp_w << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_t << endl;
    
    double stima_e = blk_av[ie]/blk_norm; // Total Energy
    glob_av[ie]  += stima_e;
    glob_av2[ie] += stima_e * stima_e;
    err_e = Statistical_uncertainty(glob_av[ie]/(double)iblk,glob_av2[ie]/(double)iblk,iblk);
    Etot_w << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_e << endl;
    
    double stima_p = blk_av[ip]/blk_norm; // Pressure
    glob_av[ip]  += stima_p;
    glob_av2[ip] += stima_p * stima_p;
    err_p = Statistical_uncertainty(glob_av[ip]/(double)iblk,glob_av2[ip]/(double)iblk,iblk);
    Pres_w << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_p << endl;
    
    Epot_w.close();
    Ekin_w.close();
    Temp_w.close();
    Etot_w.close();
    Pres_w.close();
     */
}


void ConfFinal(void){ //Write final configuration
    ofstream WriteConf;
    
    //if (readOld) {
        cout << "Print old final configuration to file old.final " << endl << endl;
        WriteConf.open("old.final");

        for (int i = 0; i < npart; ++i) {
            WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
        }
        WriteConf.close();
    //}

    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");

    for (int i=0; i<npart; ++i){
        WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    }
    WriteConf.close();
    return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r) {  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Verlet_step (double r, double old, double force, double dt) {
    return Pbc( 2.0 * r - old + force * dt*dt );
}

double Statistical_uncertainty (double average, double squared, double total_number) {
    if (total_number < 2)
        return 0;
    else
        return sqrt((squared - average*average)/(total_number - 1));
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
