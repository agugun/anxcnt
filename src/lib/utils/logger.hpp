/**
 * @file logger.hpp
 * @brief Consolidated Telemetry & Observation Layer.
 */
#pragma once
#include "lib/interfaces.hpp"
#include "lib/utils/config_reader.hpp"
#include "lib/utils/io.hpp"
#include "lib/utils/path.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <functional>
#include <memory>

namespace utl {

class StandardLogger : public top::IObserver {
public:
    using FieldExtractor = std::function<std::vector<double>(const top::IState&)>;

    struct Config {
        std::string export_dir = "exports";
        std::string project_name = "sim";
        bool enable_vtk = false;
        bool enable_csv = false;
        int terminal_freq = 10;
        int export_freq = 50;
    };

    StandardLogger(const ConfigReader& config) {
        params.export_dir = config.get("export_dir", "exports");
        params.project_name = config.get("project_name", "sim");
        params.enable_vtk = config.get("enable_vtk", 0) != 0;
        params.enable_csv = config.get("enable_csv", 0) != 0;
        params.terminal_freq = config.get("log_frequency", 10);
        params.export_freq = config.get("export_frequency", 50);
    }

    void set_grid(int nx, int ny = 1, int nz = 1, double dx = 1.0, double dy = 1.0, double dz = 1.0) {
        this->nx = nx; this->ny = ny; this->nz = nz;
        this->dx = dx; this->dy = dy; this->dz = dz;
    }

    void add_field(const std::string& name, FieldExtractor extractor) {
        fields.push_back({name, extractor});
    }

    // IObserver implementation
    void on_simulation_start(const top::IGrid& grid) override {
        std::cout << "--- Simulation START ---\n";
        step = 0;
    }

    void on_step_complete(double t, int step_idx, const top::IState& state) override {
        step = step_idx;
        if (step % params.terminal_freq == 0) {
            log_terminal(t, state);
        }
        if (step % params.export_freq == 0) {
            export_data(t, state);
        }
    }

    void on_simulation_end() override {
        std::cout << "--- Simulation END ---\n";
    }

private:
    void log_terminal(double t, const top::IState& state) {
        std::cout << "Time: " << std::fixed << std::setprecision(4) << t 
                  << " | Step: " << step << "\n";
    }

    void export_data(double t, const top::IState& state) {
        if (!params.enable_vtk && !params.enable_csv) return;

        std::vector<VTKField> vtk_fields;
        for (const auto& f : fields) {
            vtk_fields.push_back({f.name, f.extractor(state)});
        }

        if (params.enable_vtk) {
            char fname[100];
            std::sprintf(fname, "step_%04d.vti", step);
            std::string path = PathUtility::prepare_export_path(params.export_dir, params.project_name, "vtk", fname);
            
            if (nz > 1) VTKExporter::export_vti_multi_3d(path, nx, ny, nz, dx, dy, dz, vtk_fields);
            else if (ny > 1) VTKExporter::export_vti_multi_2d(path, nx, ny, dx, dy, vtk_fields);
            else VTKExporter::export_vti_multi_1d(path, nx, dx, vtk_fields);
        }

        if (params.enable_csv) {
            char fname[100];
            std::sprintf(fname, "step_%04d.csv", step);
            std::string path = PathUtility::prepare_export_path(params.export_dir, params.project_name, "csv", fname);
            CSVExporter::export_snapshot(path, nx, ny, nz, dx, dy, dz, vtk_fields);
        }
    }

    struct FieldEntry {
        std::string name;
        FieldExtractor extractor;
    };

    Config params;
    std::vector<FieldEntry> fields;
    int nx=1, ny=1, nz=1;
    double dx=1.0, dy=1.0, dz=1.0;
    int step = 0;
};

} // namespace utl
