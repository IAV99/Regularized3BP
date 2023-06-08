#include "Header.h"
#include "Integrator.h"

void evolution(const parameters& prm, const double t, const Array& x, Array* dx) {
    double mu = prm.mu;
    double mu_1 = prm.mu_1;
    double mu_2 = 1 - mu_1;
    double E = prm.E;

    double x1 = x[0], x2 = x[1], y1 = x[2], y2 = x[3], p1 = x[4], p2 = x[5], q1 = x[6], q2 = x[7];

    triple_ ai = ai_(x1, x2, y1, y2);
    triple_ bi = bi_(x1, x2, y1, y2);
    triple_ ci = ci_(x1, x2, y1, y2, p1, p2, q1, q2);
    triple_ di = di_(x1, x2, y1, y2, p1, p2, q1, q2);
    quad_ ki = ki_(ai, bi);

    double a0 = ai.i0, a1 = ai.i1, a2 = ai.i2, b0 = bi.i0, b1 = bi.i1, b2 = bi.i2;
    double c0 = ci.i0, c1 = ci.i1, c2 = ci.i2, d0 = di.i0, d1 = di.i1, d2 = di.i2;
    double k0 = ki.i0, k1 = ki.i1, k2 = ki.i2, k3 = ki.i3;

    double ro0 = a2 * a2 + b2 * b2;
    double ro1 = a1 * a1 + b1 * b1;
    double ro2 = a0 * a0 + b0 * b0;

    double L0 = c0 * c0 + d0 * d0;
    double L1 = c1 * c1 + d1 * d1;
    double L2 = c2 * c2 + d2 * d2;

    double dro0_dx1 = 2 * (y1 * a2 + y2 * b2);
    double dro0_dx2 = 2 * (-y2 * a2 + y1 * b2);
    double dro0_dy1 = 2 * (x1 * a2 + x2 * b2);
    double dro0_dy2 = 2 * (-x2 * a2 + x1 * b2);

    double dro1_dx1 = 2 * (x1 * a1 + x2 * b1);
    double dro1_dx2 = 2 * (-x2 * a1 + x1 * b1);
    double dro1_dy1 = 2 * (y1 * a1 + y2 * b1);
    double dro1_dy2 = 2 * (-y2 * a1 + y1 * b1);

    double dro2_dx1 = 2 * (x1 * a0 + x2 * b0);
    double dro2_dx2 = 2 * (-x2 * a0 + x1 * b0);
    double dro2_dy1 = -2 * (y1 * a0 + y2 * b0);
    double dro2_dy2 = 2 * (y2 * a0 - y1 * b0);

    double dL0_dx1 = -2 * (q1 * c0 + q2 * d0);
    double dL0_dx2 = -2 * (q2 * c0 - q1 * d0);
    double dL0_dy1 = 2 * (p1 * c0 + p2 * d0);
    double dL0_dy2 = 2 * (p2 * c0 - p1 * d0);

    double dL0_dp1 = 2 * (y1 * c0 - y2 * d0);
    double dL0_dp2 = 2 * (y2 * c0 + y1 * d0);
    double dL0_dq1 = 2 * (-x1 * c0 + x2 * d0);
    double dL0_dq2 = -2 * (x2 * c0 + x1 * d0);

    double dL1_dx1 = 2 * (q1 * c1 + q2 * d1);
    double dL1_dx2 = 2 * (q2 * c1 - q1 * d1);
    double dL1_dy1 = 2 * (p1 * c1 + p2 * d1);
    double dL1_dy2 = -2 * (-p2 * c1 + p1 * d1);

    double dL1_dp1 = 2 * (y1 * c1 - y2 * d1);
    double dL1_dp2 = 2 * (y2 * c1 + y1 * d1);
    double dL1_dq1 = 2 * (x1 * c1 - x2 * d1);
    double dL1_dq2 = 2 * (x2 * c1 + x1 * d1);

    double dL2_dx1 = 2 * (p1 * c2 + p2 * d2);
    double dL2_dx2 = 2 * (p2 * c2 - p1 * d2);
    double dL2_dy1 = -2 * (q1 * c2 + q2 * d2);
    double dL2_dy2 = 2 * (-q2 * c2 + q1 * d2);

    double dL2_dp1 = 2 * (x1 * c2 - x2 * d2);
    double dL2_dp2 = 2 * (x2 * c2 + x1 * d2);
    double dL2_dq1 = 2 * (-y1 * c2 + y2 * d2);
    double dL2_dq2 = -2 * (y2 * c2 + y1 * d2);

    double ero1 = L0 / mu_1 / 32 - ro0 * (1 - mu) * mu_1 - mu * mu_1 * mu_2 * ro2 - E * ro2 * ro0;
    double ero2 = L1 / mu_2 / 32 - ro0 * (1 - mu) * mu_2 - mu * mu_1 * mu_2 * ro1 - E * ro1 * ro0;
    double ero0 = mu * L2 / (1 - mu) / 32 - (1 - mu) * (ro1 * mu_1 + ro2 * mu_2) - E * ro1 * ro2;

    (*dx)[0] = (ro0 * mu * dL2_dp1 / (1 - mu) + ro1 * dL0_dp1 / mu_1 + ro2 * dL1_dp1 / mu_2) / 32;
    (*dx)[1] = (ro0 * mu * dL2_dp2 / (1 - mu) + ro1 * dL0_dp2 / mu_1 + ro2 * dL1_dp2 / mu_2) / 32;
    (*dx)[2] = (ro0 * mu * dL2_dq1 / (1 - mu) + ro1 * dL0_dq1 / mu_1 + ro2 * dL1_dq1 / mu_2) / 32;
    (*dx)[3] = (ro0 * mu * dL2_dq2 / (1 - mu) + ro1 * dL0_dq2 / mu_1 + ro2 * dL1_dq2 / mu_2) / 32;
    (*dx)[4] = -((ro0 * mu * dL2_dx1 / (1 - mu) + ro1 * dL0_dx1 / mu_1 + ro2 * dL1_dx1 / mu_2) / 32 + ero1 * dro1_dx1 + ero2 * dro2_dx1 + ero0 * dro0_dx1);
    (*dx)[5] = -((ro0 * mu * dL2_dx2 / (1 - mu) + ro1 * dL0_dx2 / mu_1 + ro2 * dL1_dx2 / mu_2) / 32 + ero1 * dro1_dx2 + ero2 * dro2_dx2 + ero0 * dro0_dx2);
    (*dx)[6] = -((ro0 * mu * dL2_dy1 / (1 - mu) + ro1 * dL0_dy1 / mu_1 + ro2 * dL1_dy1 / mu_2) / 32 + ero1 * dro1_dy1 + ero2 * dro2_dy1 + ero0 * dro0_dy1);
    (*dx)[7] = -((ro0 * mu * dL2_dy2 / (1 - mu) + ro1 * dL0_dy2 / mu_1 + ro2 * dL1_dy2 / mu_2) / 32 + ero1 * dro1_dy2 + ero2 * dro2_dy2 + ero0 * dro0_dy2);
    (*dx)[8] = ro0 * ro1 * ro2;
    //cout << K_xypq(L0, L1, L2, ro0, ro1, ro2, prm) << endl;
}

int main() {
    // parameters
    double mu = 0.0001;
    double e_1 = 0.7, omega_d = 0, a = 0.3, mu_1 = 0.5, sigma = 0.8, e_2 = e2(a, mu_1, sigma, e_1);

    cout << "use new parameters? [y/n]" << endl;
    char answ;
    cin >> answ;
    if (answ == 'y') {
        cout << "input mu, a, sigma, mu_1" << endl;
        cin >> mu >> a >> sigma >> mu_1;
        cout << "input omega_d" << endl;
        cin >> omega_d;
        cout << "[ e1_min , e1_max ] = [ " << e1_min(a, mu_1, sigma) << " , " << e1_max(a, mu_1, sigma) << " ]" << endl;
        cout << "input e1" << endl;
        cin >> e_1;
        while (e_1 > e1_max(a, mu_1, sigma) || e_1 < e1_min(a, mu_1, sigma)) {
            cout << "e1 not in [ e1_min , e1_max ] = [ " << e1_min(a, mu_1, sigma) << " , " << e1_max(a, mu_1, sigma) << " ]" << endl;
            cout << "input e1" << endl;
            cin >> e_1;
        }
        e_2 = e2(a, mu_1, sigma, e_1);

        ofstream out;
        out.open("parameters.txt");
        out << mu << " " << a << " " << sigma << " " << mu_1 << endl;
        out << omega_d << " " << e_1 << endl;
        out.close();

        while ((a * (1 - e_1 * e_1) + 1 / a * (1 - e_2 * e_2)) - 2 * (1 - e_1 * e_2 * cos(omega_d)) < 0) {
            cout << "parameters: a = " << a << "; mu_1 = " << mu_1 << "; sigma = " << sigma << endl;
            cout << "orbits intersect!" << endl;
            cout << "input omega_d" << endl;
            cin >> omega_d;
            cout << "[ e1_min , e1_max ] = [ " << e1_min(a, mu_1, sigma) << " , " << e1_max(a, mu_1, sigma) << " ]" << endl;
            cout << "input e1" << endl;
            cin >> e_1;
            while (e_1 > e1_max(a, mu_1, sigma) || e_1 < e1_min(a, mu_1, sigma)) {
                cout << "e1 not in [ e1_min , e1_max ] = [ " << e1_min(a, mu_1, sigma) << " , " << e1_max(a, mu_1, sigma) << " ]" << endl;
                cout << "input e1" << endl;
                cin >> e_1;
            }
            e_2 = e2(a, mu_1, sigma, e_1);
        }
    }
    else {
        ifstream in("parameters.txt");
        if (in.is_open()) {
            in >> mu >> a >> sigma >> mu_1;
            in >> omega_d >> e_1;
            in.close();
        }

    }

    double mu_2 = (double)1 - mu_1;
    cout << "parameters: a = " << a << "; mu_1 = " << mu_1 << "; sigma = " << sigma << endl;
    cout << "start from: e1 = " << e_1 << "; e2 = " << e_2 << "; omega_d = " << omega_d << endl;

    Integrator integrator;
    Array x;
    //initial condition vector

    vect_2d X0 = r1_0(a, e_1, omega_d), Y0 = r2_0(1, e_2), P0 = p1_0(a, mu_1, e_1, omega_d), Q0 = p2_0(1, mu_2, e_2, omega_d); 
    complex_ x0 = x0_(X0, Y0), y0 = y0_(X0, Y0), p0 = p0_(x0, y0, P0, Q0), q0 = q0_(x0, y0, P0, Q0);

    x[0] = x0.r;
    x[1] = x0.im;
    x[2] = y0.r;
    x[3] = y0.im;
    x[4] = p0.r;
    x[5] = p0.im;
    x[6] = q0.r;
    x[7] = q0.im;
    x[8] = 0;

    double x1 = x[0], x2 = x[1], y1 = x[2], y2 = x[3], p1 = x[4], p2 = x[5], q1 = x[6], q2 = x[7];

    triple_ ai = ai_(x1, x2, y1, y2);
    triple_ bi = bi_(x1, x2, y1, y2);
    triple_ ci = ci_(x1, x2, y1, y2, p1, p2, q1, q2);
    triple_ di = di_(x1, x2, y1, y2, p1, p2, q1, q2);
    quad_ ki = ki_(ai, bi);

    double a0 = ai.i0, a1 = ai.i1, a2 = ai.i2, b0 = bi.i0, b1 = bi.i1, b2 = bi.i2;
    double c0 = ci.i0, c1 = ci.i1, c2 = ci.i2, d0 = di.i0, d1 = di.i1, d2 = di.i2;
    double k0 = ki.i0, k1 = ki.i1, k2 = ki.i2, k3 = ki.i3;

    double ro00 = a2 * a2 + b2 * b2;
    double ro10 = a1 * a1 + b1 * b1;
    double ro20 = a0 * a0 + b0 * b0;

    double L00 = c0 * c0 + d0 * d0;
    double L10 = c1 * c1 + d1 * d1;
    double L20 = c2 * c2 + d2 * d2;

    double E = K_xypq(L00, L10, L20, ro00, ro10, ro20, { mu, a, mu_1, sigma, 0 }) / (ro00 * ro10 * ro20);
    //double E = ((L00 * ro10 / mu_1 + L10 * ro20 / mu_2 + L20 * ro00 * mu / (1 - mu)) / 32 - (1 - mu) * ro00 * (mu_1 * ro10 + mu_2 * ro20) - mu * mu_1 * mu_2 * ro10 * ro20) / (ro10 * ro20 * ro00);
    //double E = E0_(P0, Q0, ro00, ro10, ro20, mu, mu_1);

    double H0 = H_pr(P0, Q0 , X0, Y0, { mu, a, mu_1, sigma, 0 });
    double G0 = G_pr(P0, Q0, X0, Y0);

    //double E = H0;
    const parameters p = { mu, a, mu_1, sigma, E };

    //integrator parameters
    double t = 0;
    double h = 0.00001;

    double t_end = 100000 / mu;
    double t_save = 1;
    int save_counter = 1000000, counter = 0;

    ofstream out;
    out.open("results_for_pr.txt");
    out << setprecision(12) << p.mu << ", " << p.a << ", " << p.mu_1 << ", " << p.sigma << ", " << p.E << endl;

    cout << "t,   H - H0,   G - G0,   a1,   a2,   e1,   e2,   d_omega" << endl;

    while (t < t_end) {
        integrator.rk78(&t, &x, &h, 1e-9, 1e-6, 1e-2, [&p](const double t, const Array& x, Array* dx) { evolution(p, t, x, dx); });

        vect_2d r1 = X_({ x[0], x[1] }, { x[2], x[3] });
        vect_2d r2 = Y_({ x[0], x[1] }, { x[2], x[3] });
        vect_2d p1 = P_({ x[0], x[1] }, { x[2], x[3] }, { x[4], x[5] }, { x[6], x[7] });
        vect_2d p2 = Q_({ x[0], x[1] }, { x[2], x[3] }, { x[4], x[5] }, { x[6], x[7] });
        // out -> t, H, G, a1, a2, e1, e2, omega_d
        if (counter % save_counter == 0) {
            //out << setprecision(12) << H_pr(p1, p2, r1, r2, p) << endl;

            out << setprecision(12) << x[8] << ", " << H_pr(p1, p2, r1, r2, p) << ", " << G_pr(p1, p2, r1, r2) << ", ";
            out << a_from_pr(r1, p1, mu_1) << ", " << a_from_pr(r2, p2, mu_2) << ", ";
            out << e_from_pr(r1, p1, mu_1) << ", " << e_from_pr(r2, p2, mu_2) << ", ";
            out << omega_from_pr(r1, p1, mu_1) - omega_from_pr(r2, p2, mu_2) << ", ";
            //out << r1.x << ", " << r1.y << ", " << r2.x << ", " << r2.y;
            out << endl;
            
            //output in shell
            cout << setprecision(3) << x[8] << ", " << H_pr(p1, p2, r1, r2, p) - H0 << ", " << G_pr(p1, p2, r1, r2) - G0 << ", ";
            cout << a_from_pr(r1, p1, mu_1) << ", " << a_from_pr(r2, p2, mu_2) << ", ";
            cout << e_from_pr(r1, p1, mu_1) << ", " << e_from_pr(r2, p2, mu_2) << ", ";
            cout << omega_from_pr(r1, p1, mu_1) - omega_from_pr(r2, p2, mu_2) << endl;
        }
        counter++;
    }
    out.close();
    return 0;
}
