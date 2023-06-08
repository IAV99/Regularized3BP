#pragma once

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

const double pi = 3.14159265358979323846;

struct vect_2d {
	double x = 0, y = 0;
};

struct vect_3d {
	double x = 0, y = 0, z = 0;
};

struct triple_ {
	double i0 = 0, i1 = 0, i2 = 0;
};

struct quad_ {
	double i0 = 0, i1 = 0, i2 = 0, i3 = 0;
};

/*
struct t6_ {
	double i0 = 0, i1 = 0, i2 = 0, i3 = 0, i4 = 0, i5 = 0;
};
*/

struct complex_ {
	double r = 0, im = 0;
};

struct parameters {
	double mu = 0.001, a = 0.5, mu_1 = 0.5, sigma = 0.5, E = 0;
};

// from params.cpp
double Sigma(double a, double mu_1, double e1, double e2);
double e1(double a, double mu_1, double sigma, double e2);
double e2(double a, double mu_1, double sigma, double e1);
double e1_min(double a, double mu_1, double sigma);
double e2_min(double a, double mu_1, double sigma);
double e1_max(double a, double mu_1, double sigma);
double e2_max(double a, double mu_1, double sigma);
//from kepler_pr_wald.cpp
double scal_2d(vect_2d r, vect_2d d);
double mod_2d(vect_2d r);
vect_2d r1_0(double a, double e, double omega_d);
vect_2d r2_0(double a, double e);
vect_2d p1_0(double a, double mu, double e, double omega_d);
vect_2d p2_0(double a, double mu, double e, double omega_d);
double a_from_pr(vect_2d r, vect_2d p, double mu);
double e_from_pr(vect_2d r, vect_2d p, double mu);
double true_anomaly_from_pr(vect_2d r, vect_2d p, double mu);
double omega_from_pr(vect_2d r, vect_2d p, double mu);
double H_pr(vect_2d p1, vect_2d p2, vect_2d r1, vect_2d r2, parameters prm);
double G_pr(vect_2d p1, vect_2d p2, vect_2d r1, vect_2d r2);

vect_2d compl2vect(complex_ a);
complex_ vect2compl(vect_2d a);
double module_(complex_ a);
double module_quad_(complex_ a);
complex_ mult_(complex_ a, complex_ b);
complex_ divd_(complex_ a, complex_ b);
complex_ sum_(complex_ a, complex_ b, double c = 1, double d = 1);
complex_ pwr_(complex_ a, double c);
complex_ transp_(complex_ a);

vect_2d X_(complex_ x, complex_ y);
vect_2d Y_(complex_ x, complex_ y);
vect_2d P_(complex_ x, complex_ y, complex_ p, complex_ q);
vect_2d Q_(complex_ x, complex_ y, complex_ p, complex_ q);
complex_ x0_(vect_2d X0, vect_2d Y0);
complex_ y0_(vect_2d X0, vect_2d Y0);
complex_ p0_(complex_ x0, complex_ y0, vect_2d P0, vect_2d Q0);
complex_ q0_(complex_ x0, complex_ y0, vect_2d P0, vect_2d Q0);

double K_xypq(double L0, double L1, double L2, double ro0, double ro1, double ro2, parameters prm);
double E0_(vect_2d P0, vect_2d Q0, double ro0, double ro1, double ro2, double mu, double mu1);

triple_ ai_(double x1, double  x2, double  y1, double y2);
triple_ bi_(double x1, double  x2, double  y1, double y2);
quad_ ki_(triple_ a, triple_ b);
triple_ ci_(double x1, double  x2, double  y1, double y2, double p1, double  p2, double  q1, double q2);
triple_ di_(double x1, double  x2, double  y1, double y2, double p1, double  p2, double  q1, double q2);


complex_ sqrtX_(complex_ x, complex_ y);
complex_ sqrtY_(complex_ x, complex_ y);