#include "Header.h"

double scal_2d(vect_2d r, vect_2d d) {
	return r.x * d.x + r.y * d.y;
}

double mod_2d(vect_2d r) {
	return r.x * r.x + r.y * r.y;
}

double sign(double a) {
	if (a > 0) return 1;
	if (a < 0) return -1;
	return 0;
}

// start values from kepler to (p,r)

vect_2d r1_0(double a, double e, double omega_d) {
	vect_2d r;
	r.x = a * (1 - e) * cos(omega_d);
	r.y = a * (1 - e) * sin(omega_d);
	return r;
}

vect_2d r2_0(double a, double e) {
	vect_2d r;
	r.x = a * (1 - e);
	r.y = 0;
	return r;
}

vect_2d p1_0(double a, double mu, double e, double omega_d) {
	vect_2d r;
	r.x = -mu * sin(omega_d) * sqrt((1 - e * e) / a) /(1 - e);
	r.y = mu * cos(omega_d) * sqrt((1 - e * e) / a) / (1 - e);
	return r;
}

vect_2d p2_0(double a, double mu, double e, double omega_d) {
	vect_2d r;
	r.x = 0;
	r.y = mu * sqrt((1 - e * e) / a) / (1 - e);
	return r;
}

//values from (p,r) to kepler

double a_from_pr(vect_2d r, vect_2d p, double mu) {
	return -1 / (mod_2d(p) / (mu * mu) - 2 / sqrt(mod_2d(r)));
}

double e_from_pr(vect_2d r, vect_2d p, double mu) {
	double a = a_from_pr(r, p, mu);
	return sqrt(1 - 1 / (a * mu * mu) * (mod_2d(p) * mod_2d(r) - scal_2d(r, p) * scal_2d(r, p)));
}

double true_anomaly_from_pr(vect_2d r, vect_2d p, double mu) {
	double a = a_from_pr(r, p, mu);
	double e = e_from_pr(r, p, mu);
	double cos = 1 / e * (a * (1 - e * e) / sqrt(mod_2d(r)) - 1);
	double sin = scal_2d(r, p) / (sqrt(mod_2d(r)) * mu * e) * sqrt(a * (1 - e * e));
	if (abs(sin) < 0.000001) {
		return sign(cos) * (asin(sin) - pi/2) + pi/2;
	}
	return acos(cos) * sign(sin);
}

double omega_from_pr(vect_2d r, vect_2d p, double mu) {
	double tr_a = true_anomaly_from_pr(r, p, mu);
	double sin = r.y / sqrt(mod_2d(r));
	double cos = r.x / sqrt(mod_2d(r));
	if (abs(sin) < 0.000001) {
		return sign(cos) * (asin(sin) - pi / 2) + pi / 2 - tr_a;
	}
	return acos(cos) * sign(sin) - tr_a;
}

double H_pr(vect_2d p1, vect_2d p2, vect_2d r1, vect_2d r2, parameters prm) {
	double mu_2 = 1 - prm.mu_1;
	double H01 = mod_2d(p1) / (2 * prm.mu_1) - prm.mu_1 / sqrt(mod_2d(r1));
	double H02 = mod_2d(p2) / (2 * mu_2) - mu_2 / sqrt(mod_2d(r2));
	vect_2d dr = { r1.x - r2.x, r1.y - r2.y };
	vect_2d sum_p = { p1.x + p2.x, p1.y + p2.y };
	double V = prm.mu_1 * mu_2 / (sqrt(mod_2d(dr)));
	V += - prm.mu_1 / sqrt(mod_2d(r1)) - mu_2 / sqrt(mod_2d(r2)) - 1 / (2 * (1 - prm.mu)) * mod_2d(sum_p);
	return H01 + H02 - prm.mu * V;
}

double G_pr(vect_2d p1, vect_2d p2, vect_2d r1, vect_2d r2) {
	return (r1.x * p1.y - r1.y * p1.x) + (r2.x * p2.y - r2.y * p2.x);
}

//values for complex

complex_ transp_(complex_ a) {
	complex_ res = { a.r, -a.im };
	return res;
}

vect_2d compl2vect(complex_ a) {
	return { a.r, a.im };
}

complex_ vect2compl(vect_2d a) {
	return { a.x, a.y };
}

double module_(complex_ a) {
	return sqrt(a.r * a.r + a.im * a.im);
}

double module_quad_(complex_ a) {
	return a.r * a.r + a.im * a.im;
}

complex_ mult_(complex_ a, complex_ b) {
	return { a.r * b.r - a.im * b.im, a.r * b.im + a.im * b.r };
}

complex_ divd_(complex_ a, complex_ b) {
	complex_ res = mult_(a, transp_(b));
	double d = module_quad_(b);
	return { res.r / d, res.im / d };
}

complex_ sum_(complex_ a, complex_ b, double c, double d) {
	return { d * a.r + c * b.r, d * a.im + c * b.im };
}

complex_ pwr_(complex_ a, double c) {
	double mod_ = module_(a);
	double sin_ = a.im / mod_, cos_ = a.r / mod_;
	double arg_ = 0;
	if (abs(sin_) < 0.000001) {
		arg_ = sign(cos_) * (asin(sin_) - pi / 2) + pi / 2;
	}
	else {
		arg_ = sign(sin_) * acos(cos_);
	}
	return 	{ pow(mod_, c) * cos(c * arg_), pow(mod_, c) * sin(c * arg_) };
}

//values for waldvogel

complex_ sqrtX_(complex_ x, complex_ y) {
	return { 0.5 * ((x.r * x.r - x.im * x.im) - (y.r * y.r - y.im * y.im)), x.r * x.im - y.r * y.im };
}

complex_ sqrtY_(complex_ x, complex_ y) {
	return { 0.5 * ((x.r * x.r - x.im * x.im) + (y.r * y.r - y.im * y.im)), x.r * x.im + y.r * y.im };
}

complex_ x0_(vect_2d X0, vect_2d Y0) {
	return pwr_(sum_(pwr_(vect2compl(Y0), 0.5), pwr_(vect2compl(X0), 0.5)), 0.5);
}

complex_ y0_(vect_2d X0, vect_2d Y0) {
	return pwr_(sum_(pwr_(vect2compl(Y0), 0.5), pwr_(vect2compl(X0), 0.5), -1), 0.5);
}


triple_ ai_(double x1, double  x2, double  y1, double y2) {
	double a0 = 0.5 * ((x1 * x1 - x2 * x2) - (y1 * y1 - y2 * y2));
	double a1 = 0.5 * ((x1 * x1 - x2 * x2) + (y1 * y1 - y2 * y2));
	double a2 = x1 * y1 - x2 * y2;
	return { a0, a1, a2 };
}

triple_ bi_(double x1, double  x2, double  y1, double y2) {
	double b0 = x1 * x2 - y1 * y2;
	double b1 = x1 * x2 + y1 * y2;
	double b2 = x1 * y2 + x2 * y1;
	return { b0, b1, b2 };
}

quad_ ki_(triple_ a, triple_ b) {
	double k0 = a.i0 * a.i2 - b.i0 * b.i2;
	double k1 = a.i0 * b.i2 + a.i2 * b.i0;
	double k2 = a.i1 * a.i2 - b.i1 * b.i2;
	double k3 = a.i1 * b.i2 + a.i2 * b.i1;
	return { k0, k1, k2, k3 };
}

triple_ ci_(double x1, double  x2, double  y1, double y2, double p1, double  p2, double  q1, double q2) {
	double c0 = p1 * y1 + p2 * y2 - q1 * x1 - q2 * x2;
	double c1 = p1 * y1 + p2 * y2 + q1 * x1 + q2 * x2;
	double c2 = p1 * x1 + p2 * x2 - q1 * y1 - q2 * y2;
	return { c0, c1, c2 };
}

triple_ di_(double x1, double  x2, double  y1, double y2, double p1, double  p2, double  q1, double q2) {
	double d0 = -p1 * y2 + p2 * y1 + q1 * x2 - q2 * x1;
	double d1 = -p1 * y2 + p2 * y1 - q1 * x2 + q2 * x1;
	double d2 = -p1 * x2 + p2 * x1 + q1 * y2 - q2 * y1;
	return { d0, d1, d2 };
}

quad_ gi_(double x1, double  x2, double  y1, double y2, double p1, double  p2, double  q1, double q2) {
	triple_ ai = ai_(x1, x2, y1, y2);
	triple_ bi = bi_(x1, x2, y1, y2);
	double a0 = ai.i0, a1 = ai.i1, a2 = ai.i2, b0 = bi.i0, b1 = bi.i1, b2 = bi.i2;

	double g0 = p1 * a0 + p2 * b0 + q1 * a1 + q2 * b1;
	double g1 = p1 * b0 - p2 * a0 + q1 * b1 - q2 * a1;
	double g2 = -p1 * a0 - p2 * b0 + q1 * a1 + q2 * b1;
	double g3 = -p1 * b0 + p2 * a0 + q1 * b1 - q2 * a1;
	return { g0, g1, g2, g3 };
}

complex_ p0_(complex_ x0, complex_ y0, vect_2d P0, vect_2d Q0) {
	double x1 = x0.r, x2 = x0.im;
	double y1 = y0.r, y2 = y0.im;
	double P1 = P0.x, P2 = P0.y;
	double Q1 = Q0.x, Q2 = Q0.y;

	triple_ ai = ai_(x1, x2, y1, y2);
	triple_ bi = bi_(x1, x2, y1, y2);
	quad_ gi = gi_(x1, x2, y1, y2, P1, P2, Q1, Q2);

	double g0 = gi.i0, g1 = gi.i1, g2 = gi.i2, g3 = gi.i3;

	double p1 = 2 * (x1 * g0 - x2 * g1);
	double p2 = -2 * (x1 * g1 + x2 * g0);
	return { p1, p2 };
}

complex_ q0_(complex_ x0, complex_ y0, vect_2d P0, vect_2d Q0) {
	double x1 = x0.r, x2 = x0.im;
	double y1 = y0.r, y2 = y0.im;
	double P1 = P0.x, P2 = P0.y;
	double Q1 = Q0.x, Q2 = Q0.y;

	triple_ ai = ai_(x1, x2, y1, y2);
	triple_ bi = bi_(x1, x2, y1, y2);
	quad_ gi = gi_(x1, x2, y1, y2, P1, P2, Q1, Q2);

	double g0 = gi.i0, g1 = gi.i1, g2 = gi.i2, g3 = gi.i3;

	double q1 = 2 * (y1 * g2 - y2 * g3);
	double q2 = -2 * (y1 * g3 + y2 * g2);
	return { q1, q2 };
}

vect_2d X_(complex_ x, complex_ y) {
	triple_ a = ai_(x.r, x.im, y.r, y.im);
	triple_ b = bi_(x.r, x.im, y.r, y.im);
	return { a.i0 * a.i0 - b.i0 * b.i0, 2 * a.i0 * b.i0 };
}

vect_2d Y_(complex_ x, complex_ y) {
	triple_ a = ai_(x.r, x.im, y.r, y.im);
	triple_ b = bi_(x.r, x.im, y.r, y.im);
	return { a.i1 * a.i1 - b.i1 * b.i1, 2 * a.i1 * b.i1 };
}

vect_2d P_(complex_ x, complex_ y, complex_ p, complex_ q) {
	double x1 = x.r, x2 = x.im;
	double y1 = y.r, y2 = y.im;
	double p1 = p.r, p2 = p.im;
	double q1 = q.r, q2 = q.im;

	triple_ ai = ai_(x1, x2, y1, y2);
	triple_ bi = bi_(x1, x2, y1, y2);
	triple_ ci = ci_(x1, x2, y1, y2, p1, p2, q1, q2);
	triple_ di = di_(x1, x2, y1, y2, p1, p2, q1, q2);
	quad_ ki = ki_(ai, bi);

	double a0 = ai.i0, a1 = ai.i1, a2 = ai.i2, b0 = bi.i0, b1 = bi.i1, b2 = bi.i2;
	double c0 = ci.i0, c1 = ci.i1, c2 = ci.i2, d0 = di.i0, d1 = di.i1, d2 = di.i2;
	double k0 = ki.i0, k1 = ki.i1, k2 = ki.i2, k3 = ki.i3;

	double P1 = 0.25 * (c0 * k0 - d0 * k1) / (k0 * k0 + k1 * k1);
	double P2 = 0.25 * (c0 * k1 + d0 * k0) / (k0 * k0 + k1 * k1);
	return { P1, P2 };
}

vect_2d Q_(complex_ x, complex_ y, complex_ p, complex_ q) {
	double x1 = x.r, x2 = x.im;
	double y1 = y.r, y2 = y.im;
	double p1 = p.r, p2 = p.im;
	double q1 = q.r, q2 = q.im;

	triple_ ai = ai_(x1, x2, y1, y2);
	triple_ bi = bi_(x1, x2, y1, y2);
	triple_ ci = ci_(x1, x2, y1, y2, p1, p2, q1, q2);
	triple_ di = di_(x1, x2, y1, y2, p1, p2, q1, q2);
	quad_ ki = ki_(ai, bi);

	double a0 = ai.i0, a1 = ai.i1, a2 = ai.i2, b0 = bi.i0, b1 = bi.i1, b2 = bi.i2;
	double c0 = ci.i0, c1 = ci.i1, c2 = ci.i2, d0 = di.i0, d1 = di.i1, d2 = di.i2;
	double k0 = ki.i0, k1 = ki.i1, k2 = ki.i2, k3 = ki.i3;

	double Q1 = 0.25 * (c1 * k2 - d1 * k3) / (k2 * k2 + k3 * k3);
	double Q2 = 0.25 * (c1 * k3 + d1 * k2) / (k2 * k2 + k3 * k3);
	return { Q1, Q2 };
}

double K_xypq(double L0, double L1, double L2, double ro0, double ro1, double ro2, parameters prm) {
	double res = (L0 * ro1 / prm.mu_1 + L1 * ro2 / (1 - prm.mu_1) + prm.mu * L2 * ro0 / (1 - prm.mu)) / 32;
	res -= (1 - prm.mu) * ro0 * (ro1 * prm.mu_1 + ro2 * (1 - prm.mu_1)) + prm.mu * ro1 * ro2 * prm.mu_1 * (1 - prm.mu_1);
	return res - ro0 * ro1 * ro2 * prm.E;
}

double E0_(vect_2d P0, vect_2d Q0, double ro0, double ro1, double ro2, double mu, double mu1) {
	double mu2 = 1 - mu1;
	double px1 = P0.x, py1 = P0.y, px2 = Q0.x, py2 = Q0.y;
	return 0.5 * ((px1 * px1 + py1 * py1) / mu1 + (px2 * px2 + py2 * py2) / mu2 + mu * ((px1 + px2) * (px1 + px2) + (py1 + py2) * (py1 + py2)) / (1 - mu)) - (1 - mu) * (mu1 / ro2 + mu2 / ro1) - mu * mu1 * mu2 / ro0;
}