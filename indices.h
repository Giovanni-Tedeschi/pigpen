#ifndef INDICES_H
#define INDICES_H

struct Indices {
    int rho;
    int vx;
    int vy;
    int vz;
    int P;
    int s11;
    int s12;
    int s22;
    int s13;
    int s23;
    int s33;
    // Add more as needed
};

extern Indices idx; // Declare global variable

#endif // INDICES_H