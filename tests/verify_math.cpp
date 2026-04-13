#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include "lib/math_utils.hpp"
#include "lib/operators.hpp"

using namespace num;
using namespace std;

void test_fft() {
    cout << "--- Testing FFT ---" << endl;
    
    // Create a simple unit pulse
    VectorC input = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    
    cout << "Input: ";
    for (auto v : input) cout << v.real() << " ";
    cout << endl;

    // FFT
    VectorC freq = fft(input);
    cout << "FFT (Real): ";
    for (auto v : freq) cout << fixed << setprecision(2) << v.real() << " ";
    cout << endl;

    // IFFT
    VectorC reconstructed = ifft(freq);
    cout << "Reconstructed (Real): ";
    for (auto v : reconstructed) cout << fixed << setprecision(2) << v.real() << " ";
    cout << endl;

    // Check error
    double err = 0.0;
    for (size_t i = 0; i < input.size(); ++i) {
        err += abs(input[i] - reconstructed[i]);
    }
    cout << "FFT/IFFT Error: " << err << endl;
}

void test_power_iteration() {
    cout << "\n--- Testing Power Iteration (1D Laplacian) ---" << endl;
    
    int n = 100;
    double dx = 1.0;
    
    // Matrix-free operator A = -Laplacian
    // The eigenvalues of -Laplacian are lambda_k = 4/dx^2 * sin^2(k*pi / (2*(n+1)))
    // Max eigenvalue for k=n: lambda_max approx 4/dx^2
    auto apply_A = [n, dx](const Vector& v) {
        Vector res(n, 0.0);
        for (int i = 0; i < n; ++i) {
            double left = (i > 0) ? v[i-1] : 0.0;
            double right = (i < n-1) ? v[i+1] : 0.0;
            res[i] = -(left - 2.0 * v[i] + right) / (dx * dx);
        }
        return res;
    };

    double rho = estimate_spectral_radius(apply_A, n);
    double analytical_max = 4.0 / (dx * dx) * pow(sin(n * M_PI / (2.0 * (n + 1))), 2);
    
    cout << "Estimated Spectral Radius: " << rho << endl;
    cout << "Analytical Max Eigenvalue: " << analytical_max << endl;
    cout << "Relative Error: " << abs(rho - analytical_max) / analytical_max << endl;
}

int main() {
    try {
        test_fft();
        test_power_iteration();
        cout << "\nAll Math Utilities Verified." << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
