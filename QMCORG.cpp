// Diffusion Monte Carlo program for the 3-D harmonic oscillator

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
using namespace std;

//#include "gsl.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#define GSL_HPP

namespace gsl {

    // default Mersenne Twister using mt19937 algorithm
    static gsl_rng *gslcpp_rng = gsl_rng_alloc(gsl_rng_default);

    static double ran_uniform()        // uniform deviate in [0,1)
    {
        return gsl_rng_uniform(gslcpp_rng);
    }

    static double ran_gaussian(double sigma=1.0)    // Gaussian deviate
    {
        return gsl_ran_gaussian(gslcpp_rng, sigma);
    }

} /* namespace gsl */

namespace gsl {
    struct GSL_Error {
        std::string name;
        GSL_Error(const char* c) : name(c) { }
        GSL_Error(std::string s) : name(s) { }
    };

    inline void error(const char* c)
    {
        std::cerr << c << std::endl;
        throw GSL_Error(c);
    }

} /* end namespace gsl */

const int DIM = 3;             // dimensionality of space

double V(double *r) {          // harmonic oscillator in DIM dimensions
    double rSqd = 0;
    for (int d = 0; d < DIM; d++)
        rSqd += r[d] * r[d];
    return 0.5 * rSqd;
}

double dt;                     // Delta_t set by user
double E_T;                    // target energy

// random walkers
int N;                         // current number of walkers
int N_T;                       // desired target number of walkers
double **r;                    // x,y,z positions of walkers
bool *alive;                   // is this walker alive?

void ensureCapacity(int index) {

    static int maxN = 0;       // remember the size of the array

    if (index < maxN)          // no need to expand array
        return;                // do nothing

    int oldMaxN = maxN;        // remember the old capacity
    if (maxN > 0)
        maxN *= 2;             // double capacity
    else
        maxN = 1;
    if (index > maxN - 1)      // if this is not sufficient
        maxN = index + 1;      // increase it so it is sufficient

    // allocate new storage
    double **rNew = new double* [maxN];
    bool *newAlive = new bool [maxN];
    for (int n = 0; n < maxN; n++) {
        rNew[n] = new double [DIM];
        if (n < oldMaxN) {     // copy old values into new arrays
            for (int d = 0; d < DIM; d++)
                rNew[n][d] = r[n][d];
            newAlive[n] = alive[n];
            delete [] r[n];    // release old memory
        }
    }
    delete [] r;               // release old memory
    r = rNew;                  // point r to the new memory
    delete [] alive;
    alive = newAlive;
}

// observables
double ESum;                   // accumulator for energy
double ESqdSum;                // accumulator for variance
double rMax = 4;               // max value of r to measure psi
const int NPSI = 100;          // number of bins for wave function
double psi[NPSI];              // wave function histogram

void zeroAccumulators() {
    ESum = ESqdSum = 0;
    for (int i = 0; i < NPSI; i++)
        psi[i] = 0;
}

void initialize() {
    N = N_T;                   // set N to target number specified by user
    for (int n = 0; n < N; n++) {
        ensureCapacity(n);
        for (int d = 0; d < DIM; d++)
            r[n][d] = gsl::ran_uniform() - 0.5;
        alive[n] = true;
    }
    zeroAccumulators();
    E_T = 0;                   // initial guess for the ground state energy
}

void oneMonteCarloStep(int n) {

    // Diffusive step
    for (int d = 0; d < DIM; d++)
        r[n][d] += gsl::ran_gaussian() * sqrt(dt);

    // Branching step
    double q = exp(- dt * (V(r[n]) - E_T));
    int survivors = int(q);
    if (q - survivors > gsl::ran_uniform())
        ++survivors;

    // append survivors-1 copies of the walker to the end of the array
    for (int i = 0; i < survivors - 1; i++) {
        ensureCapacity(N);
        for (int d = 0; d < DIM; d++)
            r[N][d] = r[n][d];
        alive[N] = true;
        ++N;
    }

    // if survivors is zero, then kill the walker
    if (survivors == 0)
        alive[n] = false;
}

void oneTimeStep() {

    // DMC step for each walker
    int N_0 = N;
    for (int n = 0; n < N_0; n++)
        oneMonteCarloStep(n);

    // remove all dead walkers from the arrays
    int newN = 0;
    for (int n = 0; n < N; n++)
    if (alive[n]) {
        if (n != newN) {
            for (int d = 0; d < DIM; d++)
                r[newN][d] = r[n][d];
            alive[newN] = true;
        }
        ++newN;
    }
    N = newN;

    // adjust E_T
    E_T += log(N_T / double(N)) / 10;

    // measure energy, wave function
    ESum += E_T;
    ESqdSum += E_T * E_T;
    for (int n = 0; n < N; n++) {
        double rSqd = 0;
        for (int d = 0; d < DIM; d++)
            rSqd = r[n][d] * r[n][d];
        int i = int(sqrt(rSqd) / rMax * NPSI);
        if (i < NPSI)
            psi[i] += 1;
    }
}

int main() {

    cout << " Diffusion Monte Carlo for the 3-D Harmonic Oscillator\n"
         << " -----------------------------------------------------\n";
    cout << " Enter desired target number of walkers: ";
    cin >> N_T;
    cout << " Enter time step dt: ";
    cin >> dt;
    cout << " Enter total number of time steps: ";
    int timeSteps;
    cin >> timeSteps;

    initialize();

    // do 20% of timeSteps as thermalization steps
    int thermSteps = int(0.2 * timeSteps);
    for (int i = 0; i < thermSteps; i++)
        oneTimeStep();

    // production steps
    zeroAccumulators();
    for (int i = 0; i < timeSteps; i++) {
        oneTimeStep();
    }

    // compute averages
    double EAve = ESum / timeSteps;
    double EVar = ESqdSum / timeSteps - EAve * EAve;
    cout << " <E> = " << EAve << " +/- " << sqrt(EVar / timeSteps) << endl;
    cout << " <E^2> - <E>^2 = " << EVar << endl;
    double psiNorm = 0, psiExactNorm = 0;
    double dr = rMax / NPSI;
    for (int i = 0; i < NPSI; i++) {
        double r = i * dr;
        psiNorm += pow(r, DIM-1) * psi[i] * psi[i];
        psiExactNorm += pow(r, DIM-1) * exp(- r * r);
    }
    psiNorm = sqrt(psiNorm);
    psiExactNorm = sqrt(psiExactNorm);
    ofstream file("psi.data");
    for (int i = 0; i < NPSI; i++) {
        double r = i * dr;
        file << r << '\t' << pow(r, DIM-1) * psi[i] / psiNorm << '\t'
             << pow(r, DIM-1) * exp(- r * r / 2) / psiExactNorm << '\n';
    }
    file.close();
}
