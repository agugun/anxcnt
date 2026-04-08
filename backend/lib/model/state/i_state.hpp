#pragma once
#include <vector>

namespace model::state {

/**
 * @brief The State Abstraction (Memory)
 */
class IState {
public:
    virtual ~IState() = default;

    /**
     * @brief Initializes the memory arrays at t=0.
     */
    virtual void set_initial_conditions(const std::vector<double>& values) = 0;

    /**
     * @brief Returns the current spatial arrays (Pressure, Temp, etc.).
     */
    virtual const std::vector<double>& get_data() const = 0;
};

} // namespace model::state
