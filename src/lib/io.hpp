#pragma once
#include <string>
#include <vector>
#include <fstream>
#include "lib/modules.hpp"

namespace top {

class VTKExporter {
public:
    /**
     * @brief Exports a 3D scalar field to a Legacy VTK file.
     */
    static void export_structured_3d(const std::string& filename, 
                                     const std::vector<double>& data,
                                     int nx, int ny, int nz,
                                     double dx, double dy, double dz,
                                     const std::string& scalar_name = "Scalar") {
        std::ofstream file(filename);
        if (!file.is_open()) return;

        file << "# vtk DataFile Version 3.0\n";
        file << "NumPhys 3D Export\n";
        file << "ASCII\n";
        file << "DATASET STRUCTURED_POINTS\n";
        file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
        file << "ORIGIN 0 0 0\n";
        file << "SPACING " << dx << " " << dy << " " << dz << "\n";
        file << "POINT_DATA " << nx * ny * nz << "\n";
        file << "SCALARS " << scalar_name << " double\n";
        file << "LOOKUP_TABLE default\n";

        for (const auto& val : data) {
            file << val << "\n";
        }

        file.close();
    }

    /**
     * @brief Exports a 2D scalar field to a Legacy VTK file (as a Z=1 plane).
     */
    static void export_structured_2d(const std::string& filename,
                                     const std::vector<double>& data,
                                     int nx, int ny,
                                     double dx, double dy,
                                     const std::string& scalar_name = "Scalar") {
        export_structured_3d(filename, data, nx, ny, 1, dx, dy, 1.0, scalar_name);
    }

    /**
     * @brief Exports a 3D scalar field to a VTK XML Image Data (.vti) file.
     */
    static void export_vti_3d(const std::string& filename, 
                             const std::vector<double>& data,
                             int nx, int ny, int nz,
                             double dx, double dy, double dz,
                             const std::string& scalar_name = "Scalar") {
        std::ofstream file(filename);
        if (!file.is_open()) return;

        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "  <ImageData WholeExtent=\"0 " << nx-1 << " 0 " << ny-1 << " 0 " << nz-1 
             << "\" Origin=\"0 0 0\" Spacing=\"" << dx << " " << dy << " " << dz << "\">\n";
        file << "    <Piece Extent=\"0 " << nx-1 << " 0 " << ny-1 << " 0 " << nz-1 << "\">\n";
        file << "      <PointData Scalars=\"" << scalar_name << "\">\n";
        file << "        <DataArray type=\"Float64\" Name=\"" << scalar_name << "\" format=\"ascii\">\n";

        for (const auto& val : data) {
            file << val << " ";
        }

        file << "\n        </DataArray>\n";
        file << "      </PointData>\n";
        file << "    </Piece>\n";
        file << "  </ImageData>\n";
        file << "</VTKFile>\n";

        file.close();
    }

    /**
     * @brief Exports a 2D scalar field to a VTK XML Image Data (.vti) file.
     */
    static void export_vti_2d(const std::string& filename,
                             const std::vector<double>& data,
                             int nx, int ny,
                             double dx, double dy,
                             const std::string& scalar_name = "Scalar") {
        export_vti_3d(filename, data, nx, ny, 1, dx, dy, 1.0, scalar_name);
    }

    /**
     * @brief Exports a 1D scalar field to a VTK XML Image Data (.vti) file.
     */
    static void export_vti_1d(const std::string& filename,
                             const std::vector<double>& data,
                             int nx, double dx,
                             const std::string& scalar_name = "Scalar") {
        export_vti_3d(filename, data, nx, 1, 1, dx, 1.0, 1.0, scalar_name);
    }
};

} // namespace top
