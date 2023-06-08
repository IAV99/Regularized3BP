#include "Header.h"

double Sigma(double a, double mu_1, double e1, double e2) {
	return (mu_1 * sqrt(a * (1 - e1 * e1)) + (1 - mu_1) * sqrt(1 - e2 * e2)) / (1 - mu_1 * (1 - sqrt(a)));
}

double e1(double a, double mu_1, double sigma, double e2) {
	double sigma_1 = mu_1 * sqrt(a) / (1 - mu_1 * (1 - sqrt(a)));
	double s = (sigma - (1 - sigma_1) * sqrt(1 - e2 * e2)) / sigma_1;
	return sqrt(1 - s * s);
}

double e2(double a, double mu_1, double sigma, double e1) {
	double sigma_1 = mu_1 * sqrt(a) / (1 - mu_1 * (1 - sqrt(a)));
	double s = (sigma - sigma_1 * sqrt(1 - e1 * e1)) / (1 - sigma_1);
	return sqrt(1 - s * s);
}

double e1_min(double a, double mu_1, double sigma) {
	double sigma_1 = mu_1 * sqrt(a) / (1 - mu_1 * (1 - sqrt(a)));
	double sigma_2 = 1 - sigma_1;
	if (sigma_1 > sigma_2) {
		if (sigma <= sigma_1) {
			return sqrt(1 - (sigma / sigma_1) * (sigma / sigma_1));
		}
		return 0;
	}
	if (sigma <= sigma_1) {
		return sqrt(1 - (sigma / sigma_1) * (sigma / sigma_1));
	}
	return 0;
}

double e2_min(double a, double mu_1, double sigma) {
	double sigma_1 = mu_1 * sqrt(a) / (1 - mu_1 * (1 - sqrt(a)));
	double sigma_2 = 1 - sigma_1;
	if (sigma_1 > sigma_2) {
		if (sigma <= sigma_2) {
			return sqrt(1 - (sigma / sigma_2) * (sigma / sigma_2));
		}
		return 0;
	}
	if (sigma <= sigma_2) {
		return sqrt(1 - (sigma / sigma_2) * (sigma / sigma_2));
	}
	return 0;
}

double e1_max(double a, double mu_1, double sigma) {
	double sigma_1 = mu_1 * sqrt(a) / (1 - mu_1 * (1 - sqrt(a)));
	double sigma_2 = 1 - sigma_1;
	if (sigma_1 > sigma_2) {
		if (sigma <= sigma_2) {
			return 1;
		}
		return sqrt(1 - ((sigma - sigma_2) / sigma_1) * ((sigma - sigma_2) / sigma_1));
	}
	if (sigma <= sigma_2) {
		return 1;
	}
	return sqrt(1 - ((sigma - sigma_2) / sigma_1) * ((sigma - sigma_2) / sigma_1));
}

double e2_max(double a, double mu_1, double sigma) {
	double sigma_1 = mu_1 * sqrt(a) / (1 - mu_1 * (1 - sqrt(a)));
	double sigma_2 = 1 - sigma_1;
	if (sigma_1 > sigma_2) {
		if (sigma <= sigma_1) {
			return 1;
		}
		return sqrt(1 - ((sigma - sigma_1) / sigma_2) * ((sigma - sigma_1) / sigma_2));
	}
	if (sigma <= sigma_1) {
		return 1;
	}
	return sqrt(1 - ((sigma - sigma_1) / sigma_2) * ((sigma - sigma_1) / sigma_2));
}