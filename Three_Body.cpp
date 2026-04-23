/////////////////////////////////////////////////////////////////////////
/////////////////// 3 Body Problem - 3D - RK4 & RKCK ////////////////////
/////////////////////////////////////////////////////////////////////////
//
// Lorenzo Monti
//
// 03/2025
//

#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <fstream>
using namespace std;

//Choose method
#define RK_4    1
#define RKCK    2

#define METHOD RKCK

#define _USE_MATH_DEFINES
#define NMAX_EQ 256                     // Maximum number of equations

// Physical constants
const double G = 6.67430e-11;           // Cost gravita. (m^3 kg^-1 s^-2)
const double M_sun = 1.989e30;          // Mass of the Sun (kg)
const double M_earth = 5.972e24;        // Mass of the Earth (kg)
const double M_jupiter = 1.898e27;      // Mass of Jupiter (kg)
const double AU = 1.496e11;             // Astronomical Unit (m)
const double YEAR = 365.25 * 24 * 3600; // Year in seconds

// Function prototypes
void dYdt(double t, double *Y, double *R);
double Distance(double x1, double y1, double z1, double x2, double y2, double z2);
void RK4(double t, double *Y, void (*dYdt)(double, double *, double *), double dt, int neq);
void RK_CK(double t, double *Y, void (*dYdt)(double, double *, double *), double h, double tol, double &h_out, int neq);
double TotalEnergy(double *Y);


// Main function
int main() {

    int neq = 18; // 18 equations (x, y, z, vx, vy, vz for each body)
    double tb, te, t;
    double dt, nsteps;
    double h, tol, h_out;

    // Initial velocities
    double vy_Earth = 2*M_PI * AU / YEAR;  
    double vy_Jupiter = 2*M_PI * 5.2 * AU / YEAR / sqrt(5.2 * 5.2 * 5.2); 

    // Initial conditions (positions in meters, velocities in m/s)
    double Y[18] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,                 // Sun
        AU, 0.0, 0.0, 0.0, vy_Earth, 0.0,             // Earth
        5.2 * AU, 0.0, 0.0, 0.0, vy_Jupiter, 0.0      // Jupiter
    };

    tb = 0.0;                   // Start time
    te = 13.0 * YEAR;           // Simulate 13 years
#if METHOD == RK_4
    nsteps = 1e4;               // Number of steps
    dt = (te - tb) / nsteps;    // Time step
#elif METHOD == RKCK
    h = 1000.0;                 // Initial step size for RKCK method
    tol = 1e-6;                 // Tolerance 
    h_out = h;                  // Output step size
#endif
    t = tb;

    ofstream fdata("3body.dat");
    ofstream fenergy("3bodyenergy.dat");
    fdata << scientific << setprecision(6);

    int last_percent = -1; // To track the last printed percentage

    while (t < te) {
        #if METHOD == RK_4
        RK4(t, Y, dYdt, dt, neq);
        t += dt;
        #elif METHOD == RKCK
        RK_CK(t, Y, dYdt, h, tol, h_out, neq);
        t += h_out;
        #endif

        // Write time and positions (x, y, z) of Sun, Earth and Jupiter
        fdata << t / YEAR << " " 
              << Y[0] / AU << " " << Y[1] / AU << " " << Y[2] / AU << " "
              << Y[6] / AU << " " << Y[7] / AU << " " << Y[8] / AU << " "
              << Y[12] / AU << " " << Y[13] / AU << " " << Y[14] / AU 
              << endl;

        // Write time and total energy of the system
        fenergy << t / YEAR << " " << TotalEnergy(Y) << endl;

        // Visualize the progress of the simulation at the terminal.
        int percent = (t / te) * 100;
        if (percent > last_percent) {
            last_percent = percent;
            cout << "\r" << percent << "% completed " << flush;}
    }

    fenergy.close();
    fdata.close();
    cout << "-> The orbital data has been written to 3body.dat" << endl;

    return 0;
}


// Function to calculate the distance
double Distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
}

// Function to compute the derivatives of the system
void dYdt(double t, double *Y, double *R)
{
    // Sun
    double xS = Y[0], yS = Y[1], zS = Y[2];
    double vxS = Y[3], vyS = Y[4], vzS = Y[5];

    // Earth
    double xT = Y[6], yT = Y[7], zT = Y[8];
    double vxT = Y[9], vyT = Y[10], vzT = Y[11];

    // Jupiter
    double xJ = Y[12], yJ = Y[13], zJ = Y[14];
    double vxJ = Y[15], vyJ = Y[16], vzJ = Y[17];

    // Compute distances
    double rST = Distance(xS, yS, zS, xT, yT, zT);
    double rSJ = Distance(xS, yS, zS, xJ, yJ, zJ);
    double rTJ = Distance(xT, yT, zT, xJ, yJ, zJ);

    // Raising to the third power without pow
    double rST3 = rST * rST * rST;
    double rSJ3 = rSJ * rSJ * rSJ;
    double rTJ3 = rTJ * rTJ * rTJ;

    // Derivatives for the Sun
    R[0] = vxS;
    R[1] = vyS;
    R[2] = vzS;
    R[3] = -G * M_earth * (xS - xT) / rST3 - G * M_jupiter * (xS - xJ) / rSJ3;
    R[4] = -G * M_earth * (yS - yT) / rST3 - G * M_jupiter * (yS - yJ) / rSJ3;
    R[5] = -G * M_earth * (zS - zT) / rST3 - G * M_jupiter * (zS - zJ) / rSJ3;

    // Derivatives for the Earth
    R[6] = vxT;
    R[7] = vyT;
    R[8] = vzT;
    R[9] = -G * M_sun * (xT - xS) / rST3 - G * M_jupiter * (xT - xJ) / rTJ3;
    R[10] = -G * M_sun * (yT - yS) / rST3 - G * M_jupiter * (yT - yJ) / rTJ3;
    R[11] = -G * M_sun * (zT - zS) / rST3 - G * M_jupiter * (zT - zJ) / rTJ3;

    // Derivatives for Jupiter
    R[12] = vxJ;
    R[13] = vyJ;
    R[14] = vzJ;
    R[15] = -G * M_sun * (xJ - xS) / rSJ3 - G * M_earth * (xJ - xT) / rTJ3;
    R[16] = -G * M_sun * (yJ - yS) / rSJ3 - G * M_earth * (yJ - yT) / rTJ3;
    R[17] = -G * M_sun * (zJ - zS) / rSJ3 - G * M_earth * (zJ - zT) / rTJ3;
}

// Runge-Kutta 4th-order method
void RK4(double t, double *Y, void (*Ydot)(double, double *, double *), double dt, int neq)
{
    int i;
    double k1[NMAX_EQ], k2[NMAX_EQ], k3[NMAX_EQ], k4[NMAX_EQ];
    double Y1[NMAX_EQ];

    Ydot(t, Y, k1);
    for (i = 0; i < neq; i++) Y1[i] = Y[i] + 0.5 * dt * k1[i];

    Ydot(t + 0.5 * dt, Y1, k2);
    for (i = 0; i < neq; i++) Y1[i] = Y[i] + 0.5 * dt * k2[i];

    Ydot(t + 0.5 * dt, Y1, k3);
    for (i = 0; i < neq; i++) Y1[i] = Y[i] + dt * k3[i];

    Ydot(t + dt, Y1, k4);
    for (i = 0; i < neq; i++) {
        Y[i] += dt * (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) / 6.0;
    }
}

// Runge-Kutta Cash-Karp Method
void RK_CK(double t, double *Y, void (*dYdt)(double, double *, double *),
           double h, double tol, double &h_out, int neq)
{
    int i;
    
    // Cash-Karp coefficients
    const double c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 3.0 / 5.0,
                 c5 = 1.0, c6 = 7.0 / 8.0;
    const double b1 = 37.0 / 378.0, b3 = 250.0 / 621.0,
                 b4 = 125.0 / 594.0, b6 = 512.0 / 1771.0;
    const double b1_star = 2825.0 / 27648.0,
                 b3_star = 18575.0 / 48384.0, 
                 b4_star = 13525.0 / 55296.0,
                 b5_star = 277.0 / 14336.0, b6_star = 1.0 / 4.0;

    // Memory allocation for intermediate steps
    double k1[NMAX_EQ], k2[NMAX_EQ], k3[NMAX_EQ], k4[NMAX_EQ], k5[NMAX_EQ], k6[NMAX_EQ];
    double Y_temp[NMAX_EQ];
    double Y4[NMAX_EQ], Y5[NMAX_EQ];

    while (true) {
        // Compute k1
        dYdt(t, Y, k1);
        for (i = 0; i < neq; ++i) Y_temp[i] = Y[i] + h * c2 * k1[i];

        // Compute k2
        dYdt(t + c2 * h, Y_temp, k2);
        for (i = 0; i < neq; ++i) 
            Y_temp[i] = Y[i] + h * (3.0 / 40.0 * k1[i] +
                                    9.0 / 40.0 * k2[i]);

        // Compute k3
        dYdt(t + c3 * h, Y_temp, k3);
        for (i = 0; i < neq; ++i)
            Y_temp[i] = Y[i] + h * (3.0 / 10.0 * k1[i] - 
                                    9.0 / 10.0 * k2[i] + 
                                    6.0 / 5.0  * k3[i]);

        // Compute k4
        dYdt(t + c4 * h, Y_temp, k4);
        for (i = 0; i < neq; ++i)
            Y_temp[i] = Y[i] + h * (-11.0 / 54.0 * k1[i] +
                                    5.0 / 2.0    * k2[i] - 
                                    70.0 / 27.0  * k3[i] +
                                    35.0 / 27.0  * k4[i]);

        // Compute k5
        dYdt(t + c5 * h, Y_temp, k5);
        for (i = 0; i < neq; ++i)
            Y_temp[i] = Y[i] + h * (1631.0 / 55296.0   * k1[i] +
                                    175.0 / 512.0      * k2[i] + 
                                    575.0 / 13824.0    * k3[i] +
                                    44275.0 / 110592.0 * k4[i] + 
                                    253.0 / 4096.0     * k5[i]);

        // Compute k6
        dYdt(t + c6 * h, Y_temp, k6);

        // Compute Y4 and Y5 (fourth- and fifth-order approximations)
        for (i = 0; i < neq; ++i) {
            Y4[i] = Y[i] + h * (b1 * k1[i] + b3 * k3[i] +
                                b4 * k4[i] + b6 * k6[i]);
            Y5[i] = Y[i] + h * (b1_star * k1[i] + b3_star * k3[i] +
                                b4_star * k4[i] + b5_star * k5[i] + 
                                b6_star * k6[i]);
        }

        // Compute error estimate
        double err = 0.0;
        for (i = 0; i < neq; ++i)
            err = max(err, abs(Y5[i] - Y4[i]));

        // Error control
        if (err < tol) {
            for (i = 0; i < neq; ++i) Y[i] = Y4[i];
            h_out = h;
            
            // Increase step size if error is much smaller than tolerance
            if (err < (tol * 1.e-2)) h *= 0.9 * pow(tol / err, 0.2);
            break;
        }

        // Adjust step size (reduce if error is too large)
        h *= 0.9 * pow(tol / err, 0.25);
    }
}


// Function to calculate the total energy
double TotalEnergy(double *Y)
{
    // Positions and velocities of the Sun, Earth, and Jupiter
    double xS = Y[0], yS = Y[1], zS = Y[2];
    double vxS = Y[3], vyS = Y[4], vzS = Y[5];

    double xT = Y[6], yT = Y[7], zT = Y[8];
    double vxT = Y[9], vyT = Y[10], vzT = Y[11];

    double xJ = Y[12], yJ = Y[13], zJ = Y[14];
    double vxJ = Y[15], vyJ = Y[16], vzJ = Y[17];

    // Kinetic energy
    double KE = 0.5 * (
        M_sun   * (vxS * vxS + vyS * vyS + vzS * vzS) +
        M_earth * (vxT * vxT + vyT * vyT + vzT * vzT) +
        M_jupiter * (vxJ * vxJ + vyJ * vyJ + vzJ * vzJ)
    );

    // Potential energy (considering interactions between all pairs)
    double rST = Distance(xS, yS, zS, xT, yT, zT); // Dist. Sun-Earth
    double rSJ = Distance(xS, yS, zS, xJ, yJ, zJ); // Dist. Sun-Jupiter
    double rTJ = Distance(xT, yT, zT, xJ, yJ, zJ); // Dist. Earth-Jupiter

    double PE = - G * (
        (M_sun * M_earth) / rST +
        (M_sun * M_jupiter) / rSJ +
        (M_earth * M_jupiter) / rTJ
    );

    return KE + PE; // Total energy (kinetic + potential)
}
