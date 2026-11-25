#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>

const double PI = std::acos(-1.0);

// Potencial de Stockmayer COMPLETO (con orientaciones angulares)
// Recibe u1 = cos(theta1), u2 = cos(theta2) y phi
double potencial_SM_4D(double r, double u1, double u2, double phi,
                       double sigma, double epsilon, double mu) {
    double sr = sigma / r;
    double sr6 = std::pow(sr, 6);
    double sr12 = sr6 * sr6;
    double lj = 4.0 * epsilon * (sr12 - sr6);

    // Si mu == 0 (Lennard-Jones puro), evitamos calcular raices y cosenos
    if (mu == 0.0) {
        return lj;
    }

    // sin(theta) = sqrt(1 - cos^2(theta))
    double sin1 = std::sqrt(1.0 - u1 * u1);
    double sin2 = std::sqrt(1.0 - u2 * u2);

    // Funcion angular g(theta1, theta2, phi)
    double g = 2.0 * u1 * u2 - sin1 * sin2 * std::cos(phi);

    double dipolo = (mu * mu * g) / (r * r * r);
    return lj - dipolo;
}

// Integrando 4D: ( exp(-U_S / T*) - 1 ) * r^2
double integrando_4D(double r, double u1, double u2, double phi,
                     double T, double sigma, double epsilon, double mu) {
    double U = potencial_SM_4D(r, u1, u2, phi, sigma, epsilon, mu);
    double argumento = -U / T;

    if (argumento > 700.0) {
        argumento = 700.0;
    }

    return std::expm1(argumento) * r * r;
}

// Calculo de B2(T) mediante Montecarlo 4D (r, u1, u2, phi)
double calcular_B_MC_4D(double T, double sigma, double epsilon, double mu,
                        std::mt19937_64& rng, double r_max = 15.0, long long M_muestras = 1000000) {
    double r0 = 0.25 * sigma;

    // 4 ruletas uniformes: radio, cos(theta1), cos(theta2) y phi
    std::uniform_real_distribution<double> dist_r(r0, r_max);
    std::uniform_real_distribution<double> dist_u(-1.0, 1.0);
    std::uniform_real_distribution<double> dist_phi(0.0, 2.0 * PI);

    double suma = 0.0;

    for (long long k = 0; k < M_muestras; ++k) {
        double r_az   = dist_r(rng);
        double u1_az  = dist_u(rng);
        double u2_az  = dist_u(rng);
        double phi_az = dist_phi(rng);

        suma += integrando_4D(r_az, u1_az, u2_az, phi_az, T, sigma, epsilon, mu);
    }

    double promedio = suma / static_cast<double>(M_muestras);

    // El volumen angular (2 * 2 * 2*PI = 8*PI) por (-1/4) da exactamente -2*PI
    double integral_mc = (r_max - r0) * promedio;
    double correccion_nucleo = (2.0 / 3.0) * PI * std::pow(r0, 3);

    return correccion_nucleo - 2.0 * PI * integral_mc;
}

int main() {
    std::mt19937_64 rng(42);
    double T_min = 0.5;
    double T_max = 10.0;
    double dT = 0.1;
    long long M = 1000000;

    // CASO i: mu = 0, sigma = 1.0, variando epsilon
    std::ofstream archivo1("mc_caso_i_epsilon.csv");
    archivo1 << "T,B_eps_0.5,B_eps_1.0,B_eps_1.5\n";
    std::vector<double> eps_vals = {0.5, 1.0, 1.5};
    for (double T = T_min; T <= T_max + 1e-9; T += dT) {
        archivo1 << std::fixed << std::setprecision(2) << T;
        for (double eps : eps_vals) {
            double B = calcular_B_MC_4D(T, 1.0, eps, 0.0, rng, 15.0, M);
            archivo1 << "," << std::setprecision(6) << B;
        }
        archivo1 << "\n";
    }
    archivo1.close();

    // CASO ii: mu = 0, epsilon = 1.0, variando sigma
    std::ofstream archivo2("mc_caso_ii_sigma.csv");
    archivo2 << "T,B_sig_0.5,B_sig_1.0,B_sig_1.5\n";
    std::vector<double> sig_vals = {0.5, 1.0, 1.5};
    for (double T = T_min; T <= T_max + 1e-9; T += dT) {
        archivo2 << std::fixed << std::setprecision(2) << T;
        for (double sig : sig_vals) {
            double B = calcular_B_MC_4D(T, sig, 1.0, 0.0, rng, 15.0, M);
            archivo2 << "," << std::setprecision(6) << B;
        }
        archivo2 << "\n";
    }
    archivo2.close();

    // CASO iii: sigma = 1.0, epsilon = 1.0, variando mu (Stockmayer 4D real)
    std::ofstream archivo3("mc_caso_iii_mu.csv");
    archivo3 << "T,B_mu_0.0,B_mu_0.2,B_mu_0.4,B_mu_0.6,B_mu_0.8,B_mu_1.0,B_mu_1.2\n";
    std::vector<double> mu_vals = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
    for (double T = T_min; T <= T_max + 1e-9; T += dT) {
        archivo3 << std::fixed << std::setprecision(2) << T;
        for (double mu : mu_vals) {
            double B = calcular_B_MC_4D(T, 1.0, 1.0, mu, rng, 15.0, M);
            archivo3 << "," << std::setprecision(6) << B;
        }
        archivo3 << "\n";
    }
    archivo3.close();

    std::cout << "Simulacion Montecarlo 4D finalizada con exito.\n";
    return 0;
}
