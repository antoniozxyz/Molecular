#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>

const double PI = std::acos(-1.0);

double potencial_SM(double r, double sigma, double epsilon, double mu) {
    double sr = sigma / r;
    double sr6 = std::pow(sr, 6);
    double sr12 = sr6 * sr6;
    double lj = 4.0 * epsilon * (sr12 - sr6);
    double dipolo = (mu * mu) / (r * r * r);
    return lj - dipolo;
}

double integrando(double r, double T, double sigma, double epsilon, double mu) {
    double U = potencial_SM(r, sigma, epsilon, mu);
    double argumento = -U / T;
    if (argumento > 700.0) {
        argumento = 700.0;
    }
    return std::expm1(argumento) * r * r;
}

double calcular_B_MC(double T, double sigma, double epsilon, double mu,
                     std::mt19937_64& rng, double r_max = 15.0, long long M_muestras = 1000000) {
    double r0 = 0.25 * sigma;
    std::uniform_real_distribution<double> dist_r(r0, r_max);
    double suma = 0.0;
    for (long long k = 0; k < M_muestras; ++k) {
        double r_aleatorio = dist_r(rng);
        suma += integrando(r_aleatorio, T, sigma, epsilon, mu);
    }
    double promedio = suma / static_cast<double>(M_muestras);
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
            double B = calcular_B_MC(T, 1.0, eps, 0.0, rng, 15.0, M);
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
            double B = calcular_B_MC(T, sig, 1.0, 0.0, rng, 15.0, M);
            archivo2 << "," << std::setprecision(6) << B;
        }
        archivo2 << "\n";
    }
    archivo2.close();

    // CASO iii: sigma = 1.0, epsilon = 1.0, variando mu
    std::ofstream archivo3("mc_caso_iii_mu.csv");
    archivo3 << "T,B_mu_0.0,B_mu_0.2,B_mu_0.4,B_mu_0.6,B_mu_0.8,B_mu_1.0,B_mu_1.2\n";
    std::vector<double> mu_vals = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
    for (double T = T_min; T <= T_max + 1e-9; T += dT) {
        archivo3 << std::fixed << std::setprecision(2) << T;
        for (double mu : mu_vals) {
            double B = calcular_B_MC(T, 1.0, 1.0, mu, rng, 15.0, M);
            archivo3 << "," << std::setprecision(6) << B;
        }
        archivo3 << "\n";
    }
    archivo3.close();

    std::cout << "Simulacion Montecarlo finalizada con exito.\n";
    return 0;
}
