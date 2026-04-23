 ///////////////////////////////////////////////////////////////////////////////
////////////////////// N Body Problem - RK4 & RKCK //////////////////
///////////////////////////////////////////////////////////////////////////////
// Lorenzo Monti
//
// 02/2025
//

#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <fstream>
using namespace std;

// Choose method
#define RK_4    1
#define RKCK    2

#define METHOD RKCK

#define _USE_MATH_DEFINES
#define NMAX_EQ 48  // Massimo numero di equazioni (6 per ogni corpo)

// Costanti fisiche
const double G = 6.67430e-11;           // Costante gravitazionale (m^3 kg^-1 s^-2)
const double c = 299792458.0;           // Velocita' della luce in m/s
const double M_sun = 1.989e30;          // Massa del Sole (kg)
const double M_earth = 5.972e24;        // Massa della Terra (kg)
const double M_mars = 0.64171e24;       // Massa di Marte (kg)
const double M_jupiter = 1.898e27;      // Massa di Giove (kg)
const double M_saturn = 5.683e26;       // Massa di Saturno (kg)
const double M_uranus = 8.681e25;       // Massa di Urano (kg)
const double M_neptune = 1.024e26;      // Massa di Nettuno (kg)
const double M_pluto = 1.303e22;        // Massa di Plutone (kg)
const double AU = 1.496e11;             // Unita' astronomica (m)
const double YEAR = 365.25 * 24 * 3600; // Un anno in secondi

// Prototipi delle funzioni
void dYdt(double t, double *Y, double *R);
void RK4(double t, double *Y, void (*dYdt)(double, double *, double *), double dt, int neq);
void RK_CK(double t, double *Y, void (*dYdt)(double, double *, double *), double h, double tol, double &h_out, int neq);
double Distance(double x1, double y1, double z1, double x2, double y2, double z2);
double TotalEnergy(double *Y);


// Funzione principale
int main() {

    int neq = NMAX_EQ; // 48 equazioni (x, y, z, vx, vy, vz per ciascun corpo)
    double tb, te, t;
    double dt, nsteps;
    double h, tol, h_out;

    // Velocita' iniziali
    double vy_Earth = 2 * M_PI * AU / YEAR;  
    double vy_Mars = 2 * M_PI * 1.52 * AU / YEAR / sqrt(1.52 * 1.52 * 1.52); 
    double vy_Jupiter = 2 * M_PI * 5.2 * AU / YEAR / sqrt(5.2 * 5.2 * 5.2);
    double vy_Saturn = 2 * M_PI * 9.58 * AU / YEAR / sqrt(9.58 * 9.58 * 9.58);
    double vy_Uranus = 2 * M_PI * 19.22 * AU / YEAR / sqrt(19.22 * 19.22 * 19.22);
    double vy_Neptune = 2 * M_PI * 30.05 * AU / YEAR / sqrt(30.05 * 30.05 * 30.05);

    double inc_pluto = 17.16 * M_PI / 180.0;    // Inclinazione orbitale in radianti
    double a_pluto = 39.5 * AU;                 // Semiasse maggiore
    double e_pluto = 0.2488;                    // Eccentricita' orbitale
    double r_perihelion = a_pluto * (1 - e_pluto); // Distanza al perielio
    double v_perihelion = sqrt(G * M_sun * (2.0 / r_perihelion - 1.0 / a_pluto)); 
    double vy_Pluto = v_perihelion * cos(inc_pluto);
    double vz_Pluto = v_perihelion * sin(inc_pluto);
    
    double Y[48] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,             // Sole
        AU, 0.0, 0.0, 0.0, vy_Earth, 0.0,         // Terra
        1.52 * AU, 0.0, 0.0, 0.0, vy_Mars, 0.0,   // Marte
        5.2 * AU, 0.0, 0.0, 0.0, vy_Jupiter, 0.0, // Giove
        9.58 * AU, 0.0, 0.0, 0.0, vy_Saturn, 0.0, // Saturno
        19.22 * AU, 0.0, 0.0, 0.0, vy_Uranus, 0.0,// Urano
        30.05 * AU, 0.0, 0.0, 0.0, vy_Neptune, 0.0, // Nettuno
        r_perihelion * cos(inc_pluto), 0.0, r_perihelion * sin(inc_pluto),  // Plutone
        0.0, vy_Pluto, vz_Pluto  // Velocita' iniziali
    };    

    // Parametri iniziali
    tb = 0.0;                   // Tempo di inizio
    te = 25.0 * YEAR;           // Simulazione per x anni
#if METHOD == RK_4
    nsteps = 1e4;               // Numero di passi
    dt = (te - tb) / nsteps;    // Passo di integrazione
#elif METHOD == RKCK
    h = 1000.0;
    tol = 1e-10;
    h_out = h;
#endif
    t = tb;

    ofstream fdata("nbody.dat");
    ofstream ftime("nbody_time3.dat");
    ofstream fenergy("nbody_energy.dat");
    fdata << scientific << setprecision(6);

    int last_percent = -1; // Per monitorare il progresso

    double E0 = TotalEnergy(Y); // Compute initial total energy

    while (t < te) {
        #if METHOD == RK_4
        RK4(t, Y, dYdt, dt, neq);
        t += dt;
        #elif METHOD == RKCK
        RK_CK(t, Y, dYdt, h, tol, h_out, neq);

        // --- STAMPA DEL TIME STEP ADATTO AL TOL ---
        ftime << t / YEAR << " "    // tempo
              << h_out << " "       // passo suggerito
              << tol << endl;       // tolleranza

        t += h_out;
        #endif

        // Scrivere il tempo e le posizioni (x, y, z) di tutti i corpi
        // Output su file aggiornato
        fdata   << t / YEAR << " " 
                << Y[0] / AU << " " << Y[1] / AU << " " << Y[2] / AU << " "
                << Y[6] / AU << " " << Y[7] / AU << " " << Y[8] / AU << " "
                << Y[12] / AU << " " << Y[13] / AU << " " << Y[14] / AU << " "
                << Y[18] / AU << " " << Y[19] / AU << " " << Y[20] / AU << " "
                << Y[24] / AU << " " << Y[25] / AU << " " << Y[26] / AU << " "
                << Y[30] / AU << " " << Y[31] / AU << " " << Y[32] / AU << " "
                << Y[36] / AU << " " << Y[37] / AU << " " << Y[38] / AU << " "
                << Y[42] / AU << " " << Y[43] / AU << " " << Y[44] / AU << endl;
        
        // Compute and write relative energy error
        double E = TotalEnergy(Y);
        double err_rel = (E - E0) / fabs(E0);
        fenergy << t / YEAR << " " << err_rel << endl;


        // Visualizzare il progresso della simulazione
        int percent = (t / te) * 100;
        if (percent > last_percent) {
            last_percent = percent;
            cout << "\r" << percent << "% completato " << flush;}
    }

    fdata.close();
    ftime.close();
    fenergy.close();
    cout << "-> I dati orbitali sono stati scritti su nbody.dat" << endl;

    return 0;
}


// Function to calculate the distance
double Distance(double x1, double y1, double z1,
                double x2, double y2, double z2)
{
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * 
                (y2 - y1) + (z2 - z1) * (z2 - z1));
}

// Funzione per calcolare le derivate
void dYdt(double t, double *Y, double *R) 
{
    double positions[8][3], velocities[8][3];
    double masses[8] = {M_sun, M_earth, M_mars, M_jupiter, M_saturn, M_uranus, M_neptune, M_pluto};

    for (int i = 0; i < 8; i++) {
        positions[i][0] = Y[i * 6];     // x
        positions[i][1] = Y[i * 6 + 1]; // y
        positions[i][2] = Y[i * 6 + 2]; // z
        velocities[i][0] = Y[i * 6 + 3]; // vx
        velocities[i][1] = Y[i * 6 + 4]; // vy
        velocities[i][2] = Y[i * 6 + 5]; // vz
    }

    for (int i = 0; i < 8; i++) {
        R[i * 6] = velocities[i][0];
        R[i * 6 + 1] = velocities[i][1];
        R[i * 6 + 2] = velocities[i][2];

        double ax = 0, ay = 0, az = 0;

        for (int j = 0; j < 8; j++) {
            if (i != j) {
                double dx = positions[j][0] - positions[i][0];
                double dy = positions[j][1] - positions[i][1];
                double dz = positions[j][2] - positions[i][2];
                double r = sqrt(dx * dx + dy * dy + dz * dz);
                double r3 = r * r * r;
                double factor = G * masses[j] / r3;

                ax += factor * dx;
                ay += factor * dy;
                az += factor * dz;
            }
        }

        R[i * 6 + 3] = ax;
        R[i * 6 + 4] = ay;
        R[i * 6 + 5] = az;
    }
}




// Runge-Kutta Cash-Karp Method
void RK_CK(double t, double *Y, void (*dYdt)(double, double *, double *),
           double h, double tol, double &h_out, int neq)
{
    int i;
    // Cash-Karp coefficients
    const double c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 3.0 / 5.0, c5 = 1.0,
                 c6 = 7.0 / 8.0;
    const double b1 = 37.0 / 378.0, b3 = 250.0 / 621.0, b4 = 125.0 / 594.0,
                 b6 = 512.0 / 1771.0;
    const double b1_star = 2825.0 / 27648.0, b3_star = 18575.0 / 48384.0, 
                 b4_star = 13525.0 / 55296.0, b5_star = 277.0 / 14336.0, 
                 b6_star = 1.0 / 4.0;

    // Memory allocation for intermediate steps
    double k1[NMAX_EQ], k2[NMAX_EQ], k3[NMAX_EQ], k4[NMAX_EQ],
           k5[NMAX_EQ], k6[NMAX_EQ];
    double Y_temp[NMAX_EQ];
    double Y4[NMAX_EQ], Y5[NMAX_EQ];

    while (true) {
        // Compute k1
        dYdt(t, Y, k1);
        for (i = 0; i < neq; ++i)
            Y_temp[i] = Y[i] + h * c2 * k1[i];

        // Compute k2
        dYdt(t + c2 * h, Y_temp, k2);
        for (i = 0; i < neq; ++i)
            Y_temp[i] = Y[i] + h * (3.0 / 40.0 * k1[i] + 9.0 / 40.0 * k2[i]);

        // Compute k3
        dYdt(t + c3 * h, Y_temp, k3);
        for (i = 0; i < neq; ++i)
            Y_temp[i] = Y[i] + h * (3.0 / 10.0 * k1[i] - 9.0 / 10.0 * k2[i] + 
                                    6.0 / 5.0 * k3[i]);

        // Compute k4
        dYdt(t + c4 * h, Y_temp, k4);
        for (i = 0; i < neq; ++i)
            Y_temp[i] = Y[i] + h * (-11.0 / 54.0 * k1[i] + 5.0 / 2.0 * k2[i] - 
                                    70.0 / 27.0 * k3[i] + 35.0 / 27.0 * k4[i]);

        // Compute k5
        dYdt(t + c5 * h, Y_temp, k5);
        for (i = 0; i < neq; ++i)
            Y_temp[i] = Y[i] + h * (1631.0 / 55296.0 * k1[i] +
                                    175.0 / 512.0 * k2[i] + 
                                    575.0 / 13824.0 * k3[i] + 
                                    44275.0 / 110592.0 * k4[i] + 
                                    253.0 / 4096.0 * k5[i]);

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

        // Error control: if within tolerance, accept the step
        if (err < tol) {
            for (i = 0; i < neq; ++i)
                Y[i] = Y4[i];
            h_out = h;
            break;
        }

        // Adjust step size
        h *= 0.9 * pow(tol / err, 0.25);
    }
}


// Runge-Kutta 4th-order method
void RK4(double t, double *Y, void (*Ydot)(double, double *, double *),
         double dt, int neq)
{
    int i;
    double k1[NMAX_EQ], k2[NMAX_EQ], k3[NMAX_EQ], k4[NMAX_EQ], Y1[NMAX_EQ];

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

double TotalEnergy(double *Y) {
    const int nb = 8;
    double KE = 0.0;
    double PE = 0.0;
    double masses[nb] = { M_sun, M_earth, M_mars, M_jupiter,
                          M_saturn, M_uranus, M_neptune, M_pluto };

    // Kinetic energy
    for (int i = 0; i < nb; ++i) {
        double vx = Y[i*6 + 3], vy = Y[i*6 + 4], vz = Y[i*6 + 5];
        KE += 0.5 * masses[i] * (vx*vx + vy*vy + vz*vz);
    }

    // Potential energy (all distinct pairs)
    for (int i = 0; i < nb; ++i) {
        double xi = Y[i*6], yi = Y[i*6+1], zi = Y[i*6+2];
        for (int j = i+1; j < nb; ++j) {
            double xj = Y[j*6], yj = Y[j*6+1], zj = Y[j*6+2];
            double r = Distance(xi, yi, zi, xj, yj, zj);
            PE -= G * masses[i] * masses[j] / r;
        }
    }

    return KE + PE;
}

