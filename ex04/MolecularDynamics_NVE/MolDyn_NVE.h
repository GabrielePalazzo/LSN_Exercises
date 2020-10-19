/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=5;
int n_props;
int iv, ik, it, ie, ip;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;


const int nbins = 100;
double bin_size, stima_gdir;

const char** blockOutputFileName; //[m_props];

double blk_av[m_props]; //, blk_av2[m_props],
double blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double walker[m_props];

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double histogram[nbins];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, printStep, blockNumber, seed;
bool readOld = false;
int dataBlocking;
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
//void Measure(double);
double Force(int, int);
double Pbc(double);
double Verlet_step (double, double, double, double);
double Statistical_uncertainty (double, double, double);
void Reset(int);
void Accumulate(void);
void Averages(int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
