/**
 * @file config_reader.hpp
 * @brief Unified Configuration Reader for Simulation Parameters.
 * 
 * Objective:
 * This file provides a simple, robust Key-Value parser to externalize 
 * simulation inputs (e.g., grid size, time steps, physical properties).
 * 
 * Architectural Rationale:
 * By centralizing configuration logic, we decouple simulation parameters 
 * from the source code. This allows users to modify simulation behavior 
 * without recompiling, adhering to the principle of "Separation of Concerns."
 */
#pragma once
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace utl {

class ConfigReader {
public:
    ConfigReader() = default;

    /**
     * @brief Loads a configuration file.
     * @param filename Path to the params.txt file.
     * @return true if successful.
     */
    bool load(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "[ConfigReader] Warning: Could not open " << filename << ". Using defaults.\n";
            return false;
        }

        std::string line;
        while (std::getline(file, line)) {
            // Trim whitespace
            line.erase(0, line.find_first_not_of(" \t\r\n"));
            line.erase(line.find_last_not_of(" \t\r\n") + 1);

            // Skip comments and empty lines
            if (line.empty() || line[0] == '#' || (line.size() >= 2 && line[0] == '/' && line[1] == '/')) {
                continue;
            }

            size_t eq_pos = line.find('=');
            if (eq_pos == std::string::npos) continue;

            std::string key = line.substr(0, eq_pos);
            std::string val = line.substr(eq_pos + 1);

            // Trim key and value
            key.erase(key.find_last_not_of(" \t") + 1);
            val.erase(0, val.find_first_not_of(" \t"));
            val.erase(val.find_last_not_of(" \t\r\n") + 1);

            data[key] = val;
        }
        return true;
    }

    /**
     * @brief Gets a value as a specific type.
     * @param key The parameter name.
     * @param default_val Value to return if key is missing.
     */
    template<typename T>
    T get(const std::string& key, T default_val) const {
        auto it = data.find(key);
        if (it == data.end()) return default_val;

        std::stringstream ss(it->second);
        T res;
        if (ss >> res) return res;
        return default_val;
    }

    // Specialized string getter
    std::string get(const std::string& key, const std::string& default_val) const {
        auto it = data.find(key);
        if (it == data.end()) return default_val;
        return it->second;
    }

    // String literal specialization
    std::string get(const std::string& key, const char* default_val) const {
        return get(key, std::string(default_val));
    }

    bool has(const std::string& key) const {
        return data.find(key) != data.end();
    }

    /**
     * @brief Programmatically set or override a parameter.
     */
    template<typename T>
    void set(const std::string& key, T value) {
        std::stringstream ss;
        ss << value;
        data[key] = ss.str();
    }

private:
    std::map<std::string, std::string> data;
};

} // namespace utl
