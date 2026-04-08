#include "simulation.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

/**
 * @brief Main entry point for standalone 1D Pressure Diffusivity simulation.
 */
int main() {
    // Reservoir setup
    int nx = 50;
    double dx = 20.0;    // ft
    double k = 50.0;     // md
    double phi = 0.25;   // fraction
    double mu = 2.0;     // cp
    double ct = 1e-5;    // psi^-1
    
    numerical_methods::physics_pressure::PressureSimulation sim(nx, dx, k, phi, mu, ct);
    
    // Boundary conditions (Pressure Drawdown)
    double p_left = 1500.0;  // Well bottom-hole pressure (psi)
    double p_right = 3000.0; // Initial reservoir pressure (psi)
    sim.set_boundary_conditions(p_left, p_right);
    
    // Initial condition: Uniform reservoir pressure
    std::vector<double> ic(nx, p_right);
    sim.set_initial_condition(ic);
    
    // Prepare output file in dedicated directory
    std::ofstream csv_file("output/pressure_results.csv");
    csv_file << "time,index,x,pressure\n";

    // Step configuration
    double dt = 24.0;    // 1 day (hours)
    int steps = 100;
    
    std::cout << "Starting 1D Pressure Diffusivity Simulation (Newton-Raphson)\n";
    std::cout << "Diffusivity (eta): " << 0.0002637 * k / (phi * mu * ct) << " ft^2/hr\n";
    std::cout << "Results will be saved to pressure_results.csv\n";
    
    double current_time = 0.0;
    // Initial state
    sim.write_output(csv_file, sim.getState(), current_time);

    for (int t = 1; t <= steps; ++t) {
        sim.step(dt);
        current_time += dt;

        // Save every 5 steps
        if (t % 5 == 0) {
            sim.write_output(csv_file, sim.getState(), current_time);
        }
        
        if (t % 10 == 0) {
            const auto& p = sim.get_values();
            std::cout << "Step: " << t << " Time: " << current_time 
                      << " Mid-point: " << p[nx/2] << " psi\n";
        }
    }
    
    std::cout << "Simulation Successful. Detailed data in pressure_results.csv\n";
    return 0;
}
