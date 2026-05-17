#include "simulation.hpp"
#include "lib/utils/config_reader.hpp"

using namespace mod;
using namespace sim;
using namespace utl;
using namespace num;







/**
 * @brief Observer to write st variables to CSV.
 */
class CSVObserver : public IObserver {
    std::ofstream out;
public:
    CSVObserver(const std::string& filename) {
        out.open(filename);
        out << "Time,Position,Velocity\n";
    }
    void on_step_complete(double t, int step, const IState& base_state) override {
        const auto& st = dynamic_cast<const OscillatorState&>(base_state);
        out << t << "," << st.x << "," << st.v << "\n";
    }
};

int main(int argc, char** argv) {
    std::cout << "Starting Harmonic Oscillator Simulation..." << std::endl;

    ConfigReader config;
    if (argc > 1) config.load(argv[1]);

    // 1. Initialize Specialized Simulation
    OscillatorSimulation sim;
    sim.build(config);

    // 2. Initial State
    auto st_init = sim.create_initial_state(config);

    // 3. Observer Setup
    auto observer = std::make_shared<CSVObserver>("exports/oscillator_output.csv");
    sim.add_observer(observer);

    // 4. Run Simulation
    double t_max = config.get("t_max", 10.0);
    double dt = config.get("dt", 0.05);

    sim.run(t_max + 1e-4, dt, std::move(st_init));

    std::cout << "Simulation complete. Output saved to exports/oscillator_output.csv" << std::endl;
    return 0;
}
