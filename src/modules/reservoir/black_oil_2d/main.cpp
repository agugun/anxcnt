#include <omp.h>
#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include "state.hpp"
#include "model.hpp"
#include "lib/solvers.hpp"
#include "lib/integrators.hpp"
#include "lib/io.hpp"
#include "lib/operators.hpp"
#include "lib/config_reader.hpp"

using namespace top;
using namespace mod;
using namespace mod::reservoir;
using namespace num;

// Simple constant rate producer for Black Oil
class ReservoirWellBlackOil2D : public ISourceSink {
public:
    int i, j;
    double q_total;
    
    ReservoirWellBlackOil2D(int i_v, int j_v, double q) : i(i_v), j(j_v), q_total(q) {}

    void apply(Vector& residual, Matrix* jacobian, const top::IState& state, double dt,
               std::vector<num::SparseMatrix::Entry>* sparse_entries = nullptr) override {
        const auto& b_state = dynamic_cast<const ReservoirBlackOil2DState&>(state);
        int c = b_state.idx(i, j);
        
        // Oil-only production
        residual[3 * c + 1] += q_total; 
    }
};

int main(int argc, char** argv) {
    std::string config_file = "input/reservoir_black_oil_2d.txt";
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') {
            config_file = argv[i];
            break;
        }
    }

    ConfigReader config;
    config.load(config_file);
    int num_threads = config.get("num_threads", 1);
    omp_set_num_threads(num_threads);

    bool enable_vtk = config.get("enable_vtk", 0) != 0;
    bool enable_csv = config.get("enable_csv", 0) != 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--vtk") enable_vtk = true;
        else if (arg == "--csv") enable_csv = true;
    }

    std::cout << "--- 2D Full Black Oil Simulation (3-Phase) ---" << std::endl;
    if (enable_vtk) std::cout << "VTK Export: Enabled" << std::endl;
    if (enable_csv) std::cout << "CSV Export: Enabled" << std::endl;

    // 1. Setup Grid
    int nx = config.get("nx", 10);
    int ny = config.get("ny", 10);
    double dx = config.get("dx", 100.0);
    double dy = config.get("dy", 100.0);
    Spatial2D spatial(nx, ny, dx, dy);
    auto state = std::make_shared<ReservoirBlackOil2DState>(spatial);

    // 2. Define Wells
    std::vector<std::shared_ptr<ISourceSink>> wells;
    wells.push_back(std::make_shared<ReservoirWellBlackOil2D>(nx/2, ny/2, config.get("q_well", -200.0))); 

    // 3. Setup Model
    auto model = std::make_shared<ReservoirBlackOil2DModel>(config.get("k", 100.0), config.get("phi", 0.2), config.get("h", 50.0), wells);

    // 4. Setup Simulator
    auto solver = std::make_shared<NewtonSolver>();
    auto integrator = std::make_shared<FullyImplicitIntegrator>();
    StandardSimulator sim(model, state, solver, integrator);

    // 5. Run Simulation - Units: DAYS
    double t_end = config.get("t_end", 1.0);      
    double dt = config.get("dt", 0.01); 

    std::cout << "Starting 3-Phase Black Oil Depletion for " << t_end << " days..." << std::endl;

    // Prep CSV file if enabled
    std::ofstream csv_file;
    if (enable_csv) {
        csv_file.open("exports/black_oil_history.csv");
        csv_file << "Time,AvgPressure,AvgSw,AvgSg\n";
    }

    sim.run(t_end, dt, [&](double t, const IState& s) {
        const auto& cur_state = dynamic_cast<const ReservoirBlackOil2DState&>(s);
        
        // Report/Export every 0.1 days
        if (std::abs(fmod(t, 0.1)) < 1e-7 || std::abs(fmod(t, 0.1) - 0.1) < 1e-7) {
            double avg_p = 0.0, avg_sw = 0.0, avg_sg = 0.0;
            for(int i=0; i<nx*ny; ++i) {
                avg_p += cur_state.p(i);
                avg_sw += cur_state.sw(i);
                avg_sg += cur_state.sg(i);
            }
            avg_p /= (nx*ny); avg_sw /= (nx*ny); avg_sg /= (nx*ny);

            std::cout << "Time: " << t << " days | Avg P: " << avg_p << " psi | Avg Sw: " << avg_sw << " | Avg Sg: " << avg_sg << std::endl;

            if (enable_csv && csv_file.is_open()) {
                csv_file << t << "," << avg_p << "," << avg_sw << "," << avg_sg << "\n";
            }

            if (enable_vtk) {
                // VTK Export
                std::string filename = "exports/black_oil_" + std::to_string((int)(t*100)) + ".vti";
                
                std::vector<double> p_v(nx*ny), sw_v(nx*ny), sg_v(nx*ny);
                for(int i=0; i<nx*ny; ++i) {
                    p_v[i] = cur_state.p(i);
                    sw_v[i] = cur_state.sw(i);
                    sg_v[i] = cur_state.sg(i);
                }
                
                std::vector<VTKField> export_fields = {
                    {"Pressure", p_v},
                    {"WaterSaturation", sw_v},
                    {"GasSaturation", sg_v}
                };
                VTKExporter::export_vti_multi_2d(filename, nx, ny, dx, dy, export_fields);
            }
        }
    });

    if (enable_csv && csv_file.is_open()) {
        csv_file.close();
    }

    std::cout << "Successfully completed Black Oil Simulation." << std::endl;

    return 0;
}
