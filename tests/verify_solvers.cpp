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
    Matrix A = {
        {4.0, 1.0},
        {1.0, 3.0}
    };
    Vector b = {1.0, 2.0};
    
    auto apply_A = [&A](const Vector& v) {
        Vector res(2, 0.0);
        res[0] = A[0][0] * v[0] + A[0][1] * v[1];
        res[1] = A[1][0] * v[0] + A[1][1] * v[1];
        return res;
    };

    Vector x0 = {0.0, 0.0};
    auto start = chrono::high_resolution_clock::now();
    Vector x = solver.solve_iterative(apply_A, b, x0);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, micro> duration = end - start;

    cout << "Solution: ";
    for (auto val : x) cout << fixed << setprecision(4) << val << " ";
    cout << "\nExecution Time: " << duration.count() << " us" << endl;

    double total_err = 0;
    Vector Ax = apply_A(x);
    for (int i = 0; i < 2; ++i) total_err += abs(Ax[i] - b[i]);
    
    cout << "Residual Error: " << total_err << endl;
    assert(total_err < 1e-6);
}

class SimpleState : public IState {
public:
    Vector val;
    SimpleState(double x) : val({x}) {}
    void update(const Vector& delta) override { val[0] += delta[0]; }
    unique_ptr<IState> clone() const override { return make_unique<SimpleState>(val[0]); }
};

/**
 * @brief Simple Mock Model for Newton Testing (Solves f(x) = x^2 - 2 = 0)
 */
class SquareRootModel : public IModel {
public:
    Vector evaluate_rhs(const IState& state) const override { return {}; }
    
    Vector build_residual(const IState& state, const IState& state_old, double dt) const override {
        const auto& s = dynamic_cast<const SimpleState&>(state);
        double x = s.val[0];
        return { x * x - 2.0 };
    }
    
    Matrix build_jacobian(const IState& state, double dt) const override {
        const auto& s = dynamic_cast<const SimpleState&>(state);
        double x = s.val[0];
        return { { 2.0 * x } };
    }
};

void test_newton_solver() {
    cout << "\n--- Testing Newton Solver (Solve x^2 - 2 = 0) ---" << endl;
    SquareRootModel model;
    SimpleState state(1.0); // Guess x = 1.0

    auto start = chrono::high_resolution_clock::now();
    auto result = NewtonSolver::solve_robust(model, state, 0.0);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, micro> duration = end - start;

    cout << "Converged: " << (result.converged ? "Yes" : "No") << endl;
    cout << "Iterations: " << result.iterations << endl;
    cout << "Final x: " << fixed << setprecision(6) << state.val[0] << " (Expected: 1.414214)" << endl;
    cout << "Execution Time: " << duration.count() << " us" << endl;
    
    assert(result.converged);
    assert(abs(state.val[0] - sqrt(2.0)) < 1e-4);
}

int main() {
    try {
        test_lu_solver();
        test_cg_solver();
        test_newton_solver();
        cout << "\nAll Solver Tests Passed." << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
