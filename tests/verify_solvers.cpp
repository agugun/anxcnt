#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <chrono>
#include "lib/solvers.hpp"

using namespace num;
using namespace std;

void test_lu_solver() {
    cout << "--- Testing LU Solver ---" << endl;
    LUSolver solver;

    // System:
    // 3x + 2y - z = 1
    // 2x - 2y + 4z = -2
    // -x + 0.5y - z = 0
    Matrix A = {
        {3.0, 2.0, -1.0},
        {2.0, -2.0, 4.0},
        {-1.0, 0.5, -1.0}
    };
    Vector b = {1.0, -2.0, 0.0};

    auto start = chrono::high_resolution_clock::now();
    Vector x = solver.solve_system(A, b);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, micro> duration = end - start;

    cout << "Solution: ";
    for (auto val : x) cout << fixed << setprecision(4) << val << " ";
    cout << "\nExecution Time: " << duration.count() << " us" << endl;

    // Verification: Ax - b approx 0
    double total_err = 0;
    for (int i = 0; i < 3; ++i) {
        double row_sum = 0;
        for (int j = 0; j < 3; ++j) row_sum += A[i][j] * x[j];
        total_err += abs(row_sum - b[i]);
    }
    cout << "Residual Error: " << total_err << endl;
    assert(total_err < 1e-10);
}

void test_cg_solver() {
    cout << "\n--- Testing Conjugate Gradient Solver ---" << endl;
    ConjugateGradientSolver solver;

    // Symmetric Positive Definite system:
    // 4x + y = 1
    // x + 3y = 2
    SparseMatrix A(2, 2);
    A.triplets.push_back({0, 0, 4.0});
    A.triplets.push_back({0, 1, 1.0});
    A.triplets.push_back({1, 0, 1.0});
    A.triplets.push_back({1, 1, 3.0});
    A.compress();
    
    Vector b = {1.0, 2.0};
    Vector x0 = {0.0, 0.0};
    
    auto start = chrono::high_resolution_clock::now();
    Vector x = solver.solve_system(A, b, x0);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, micro> duration = end - start;

    cout << "Solution: ";
    for (auto val : x) cout << fixed << setprecision(4) << val << " ";
    cout << "\nExecution Time: " << duration.count() << " us" << endl;

    double total_err = 0;
    Vector Ax = A.multiply(x);
    for (int i = 0; i < 2; ++i) total_err += abs(Ax[i] - b[i]);
    
    cout << "Residual Error: " << total_err << endl;
    assert(total_err < 1e-6);
}

int main() {
    try {
        test_lu_solver();
        test_cg_solver();
        cout << "\nAll Solver Tests Passed." << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
