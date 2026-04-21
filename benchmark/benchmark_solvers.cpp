/**
 * @file benchmark_solvers.cpp
 * @brief Multithreaded performance benchmarking for linear solvers.
 */
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include "lib/modules.hpp"
#include "lib/sparse.hpp"
#include "lib/solvers.hpp"

using namespace top;
using namespace num;
using namespace std;

/**
 * @brief SolverBenchmark Class
 * 
 * Objectives:
 * 1. Generate a large 3D Laplacian system to test parallel efficiency.
 * 2. Profile solving time across different thread counts.
 * 3. Report speedup and resource utilization.
 */
class SolverBenchmark {
private:
    int nx, ny, nz;
    int n_total;
    SparseMatrix A;
    Vector b;
    Vector x0;

    struct Result {
        int threads;
        double time_ms;
        double speedup;
    };
    vector<Result> results;

public:
    SolverBenchmark(int grid_size) 
        : nx(grid_size), ny(grid_size), nz(grid_size), A(0, 0) {
        n_total = nx * ny * nz;
        A = generate_3d_laplacian();
        b = Vector(n_total, 1.0); // Constant RHS
        x0 = Vector(n_total, 0.0); // Initial guess
    }

    /**
     * @brief Generates a 7-point stencil 3D Laplacian Matrix.
     */
    SparseMatrix generate_3d_laplacian() {
        vector<SparseMatrix::Entry> entries;
        entries.reserve(n_total * 7);

        auto idx = [&](int i, int j, int k) {
            return k * (nx * ny) + j * nx + i;
        };

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int c = idx(i, j, k);
                    
                    // Diagonal
                    entries.push_back({c, c, 6.0});
                    
                    // Neighbors
                    if (i > 0) entries.push_back({c, idx(i - 1, j, k), -1.0});
                    if (i < nx - 1) entries.push_back({c, idx(i + 1, j, k), -1.0});
                    if (j > 0) entries.push_back({c, idx(i, j - 1, k), -1.0});
                    if (j < ny - 1) entries.push_back({c, idx(i, j + 1, k), -1.0});
                    if (k > 0) entries.push_back({c, idx(i, j, k - 1), -1.0});
                    if (k < nz - 1) entries.push_back({c, idx(i, j, k + 1), -1.0});
                }
            }
        }
        return SparseMatrix::from_triplets(n_total, n_total, entries);
    }

    /**
     * @brief Runs the solver for a specific thread count.
     */
    double run_single(int threads) {
        omp_set_num_threads(threads);
        
        auto start = chrono::high_resolution_clock::now();
        // Use BiCGSTAB as it's the target for multithreading
        BiCGSTABSolver::solve(A, b, 1e-6, 500, false);
        auto end = chrono::high_resolution_clock::now();
        
        chrono::duration<double, milli> duration = end - start;
        return duration.count();
    }

    /**
     * @brief Executes a gradual benchmark suite from 1 to max threads.
     */
    void run_full_suite() {
        int max_threads = omp_get_max_threads();

        cout << "Executing Gradual Solver Multithreading Benchmark..." << endl;
        cout << "Problem Size: " << n_total << " unknowns (Grid: " << nx << "x" << ny << "x" << nz << ")" << endl;
        cout << "Thread Sweep: 1 to " << max_threads << endl;
        cout << "----------------------------------------------------" << endl;

        double base_time = 0;
        for (int t = 1; t <= max_threads; ++t) {
            cout << "Running with " << t << " threads... " << flush;
            double time = run_single(t);
            if (t == 1) base_time = time;
            
            results.push_back({t, time, base_time / time});
            cout << fixed << setprecision(2) << time << " ms" << endl;
        }
    }

    /**
     * @brief Prints a formatted performance report.
     */
    void print_report() {
        cout << "\n================ Performance Report ================" << endl;
        cout << setw(10) << "Threads" << setw(15) << "Time (ms)" << setw(15) << "Speedup" << endl;
        cout << "----------------------------------------------------" << endl;
        for (const auto& res : results) {
            cout << setw(10) << res.threads 
                 << setw(15) << fixed << setprecision(2) << res.time_ms 
                 << setw(15) << fixed << setprecision(2) << res.speedup << "x" << endl;
        }
        cout << "====================================================" << endl;
    }
};

int main() {
    try {
        // Size 50 gives 125,000 unknowns, enough for scaling visibility.
        SolverBenchmark bench(50);
        bench.run_full_suite();
        bench.print_report();
    } catch (const exception& e) {
        cerr << "Benchmark failed: " << e.what() << endl;
        return 1;
    }
    return 0;
}
