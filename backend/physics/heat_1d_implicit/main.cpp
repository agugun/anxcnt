#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "simulation.hpp"

int main() {
    // Problem parameters
    int nx = 101;
    double dx = 0.01;
    double alpha = 0.05;
    double dt = 0.01;
    int steps = 100;

    // Orchestrate simulation
    numerical_methods::physics_heat::HeatSimulation sim(nx, dx, alpha);

    // Initial Condition (Gaussian Pulse)
    std::vector<double> ic(nx, 0.0);
    for (int i = 0; i < nx; ++i) {
        double x_val = i * dx;
        ic[i] = std::exp(-std::pow(x_val - 0.5, 2) / 0.01);
    }
    // Boundary Zeroing
    ic[0] = 0.0;
    ic[nx-1] = 0.0;

    sim.set_initial_condition(ic);
    sim.set_boundary_conditions(0.0, 0.0);

    // Prepare output file in dedicated directory
    std::ofstream csv_file("output/heat_results.csv");
    csv_file << "time,index,x,temperature\n";

    std::cout << "Starting Simulation: " << nx << " points, " << steps << " steps." << std::endl;
    std::cout << "Results saved to: backend/output/heat_results.csv" << std::endl;

    double current_time = 0.0;
    // Initial state
    sim.write_output(csv_file, sim.getState(), current_time);

    for (int step = 0; step < steps; ++step) {
        sim.step(dt);
        current_time += dt;
        
        // Save intermediate steps for transient visualization
        // (Logging every 5 steps for smoother animation)
        if (step % 5 == 0) {
            sim.write_output(csv_file, sim.getState(), current_time);
        }

        // Print to console every 20 steps
        if (step % 20 == 0) {
            const auto& data = sim.get_values();
            std::cout << "Step: " << step << " Time: " << current_time 
                      << " Mid-point: " << data[nx/2] << std::endl;
        }
    }

    std::cout << "Simulation Successful. Detailed data in heat_results.csv" << std::endl;
    return 0;
}
