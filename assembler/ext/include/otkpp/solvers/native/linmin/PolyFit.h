
#ifndef POLYFIT_H

int cubic_minimizer(double a, double fa, double dfa,
                    double b, double fb, double dfb,
                    int extrapolate, double *alpha_min);

int quad_minimizer1(double a, double fa, double dfa,
                    double b, double fb,
                    double *alpha_min);

int quad_minimizer2(double a, double dfa,
                    double b, double dfb,
                    double *alpha_min);

#define POLYFIT_H

#endif
