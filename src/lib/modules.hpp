/**
 * @file modules.hpp
 * @brief Forwarding header for Simulation Interfaces.
 * 
 * This file is deprecated in favor of lib/interfaces.hpp.
 * It remains to support legacy include paths during the architectural migration.
 */
#pragma once
#include "lib/interfaces.hpp"

namespace top {
    // Deprecated alias for SimulationEngine if needed elsewhere
    // class StandardSimulator; // Removed in favor of SimulationEngine
}