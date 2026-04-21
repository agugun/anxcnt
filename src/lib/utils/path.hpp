#pragma once
#include <string>
#include <filesystem>
#include <iostream>

namespace utl {

namespace fs = std::filesystem;

class PathUtility {
public:
    /**
     * @brief Prepares and returns the full export path, ensuring directories exist.
     * 
     * @param base_dir The root export directory (e.g., "exports")
     * @param project_name The name of the project/case (e.g., "alpha_run")
     * @param filename The final filename (e.g., "sim_results.vti")
     * @return std::string The full path for the file.
     */
    static std::string prepare_export_path(const std::string& base_dir, 
                                           const std::string& project_name, 
                                           const std::string& format,
                                           const std::string& filename) {
        
        fs::path base(base_dir);
        fs::path project(project_name);
        fs::path fmt(format);
        fs::path dir = base / project / fmt;

        try {
            if (!fs::exists(dir)) {
                fs::create_directories(dir);
            }
        } catch (const fs::filesystem_error& e) {
            std::cerr << "[PathUtility] Error creating directory: " << e.what() << "\n";
            return filename; // Fallback to current dir if error
        }

        return (dir / filename).string();
    }

    /**
     * @brief Shortcut to get project-specific path from config.
     */
    static std::string get_path(const class ConfigReader& config, const std::string& default_filename) {
        // We'll define this after forward-declaring ConfigReader if needed, 
        // but for now let's keep it simple.
        return "";
    }
};

} // namespace utl
