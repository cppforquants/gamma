#include <iostream>
#include <cmath>

// Standard normal probability density function
double norm_pdf(double x) {
    return (1.0 / std::sqrt(2 * M_PI)) * std::exp(-0.5 * x * x);
}

// Compute d1 used in the Black-Scholes formula
double compute_d1(double S, double K, double r, double sigma, double T) {
    return (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
}

// Compute Gamma using the closed-form Black-Scholes formula
double gamma(double S, double K, double r, double sigma, double T) {
    double d1 = compute_d1(S, K, r, sigma, T);
    return norm_pdf(d1) / (S * sigma * std::sqrt(T));
}

int main() {
    double S = 100.0;     // Spot price
    double K = 100.0;     // Strike price
    double r = 0.05;      // Risk-free interest rate
    double sigma = 0.2;   // Volatility
    double T = 1.0;       // Time to maturity in years

    double gamma_val = gamma(S, K, r, sigma, T);
    std::cout << "Gamma: " << gamma_val << std::endl;

    return 0;
}
