#include "lib/modules.hpp"
#include "modules/heat_1d_explicit/state.hpp"
#include "lib/operators.hpp"

namespace mod {
using namespace top;
namespace physics_heat {

class HeatModel : public IModel {
private:
    double alpha;
    double dx;

public:
    HeatModel(double alpha, double dx) : alpha(alpha), dx(dx) {}

    Vector evaluate_rhs(const IState& state) const override {
        const auto& heat_state = dynamic_cast<const HeatState&>(state);
        Vector rhs = mop::laplace_1d(heat_state.temperatures, dx);
        
        // Multiply by alpha
        for (auto& val : rhs) val *= alpha;
        
        return rhs;
    }


    // Dummy implementations for implicit methods
    Vector build_residual(const IState&, const IState&, double) const override { return {}; }
    Matrix build_jacobian(const IState&, double) const override { return {}; }
};

} // namespace physics_heat
} // namespace mod