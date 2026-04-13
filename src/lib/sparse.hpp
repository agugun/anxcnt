/**
 * @file sparse.hpp
 * @brief Compressed Sparse Row (CSR) Matrix Infrastructure.
 * 
 * Objective:
 * Provides efficient storage and arithmetic for sparse systems, which are
 * characteristic of discretized reservoir simulations where most cells 
 * only interact with their immediate neighbors.
 * 
 * Strategic Importance:
 * Transitioning from O(N^2) dense storage to O(N) sparse storage is the 
 * single most impactful optimization for large-scale 2D and 3D simulations.
 */
#pragma once
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace num {

/**
 * @brief Compressed Sparse Row (CSR) Matrix.
 */
class SparseMatrix {
public:
    int rows, cols;
    std::vector<double> values;     // Non-zero values
    std::vector<int> col_indices;   // Column index of each non-zero
    std::vector<int> row_ptr;       // Index in 'values' where each row starts

    SparseMatrix(int r, int c) : rows(r), cols(c), row_ptr(r + 1, 0) {}

    /**
     * @brief Build CSR from a list of triplets (row, col, value).
     * Automatically handles duplicates by summing them.
     */
    struct Entry {
        int r, c;
        double v;
        bool operator<(const Entry& other) const {
            if (r != other.r) return r < other.r;
            return c < other.c;
        }
    };

    static SparseMatrix from_triplets(int r, int c, std::vector<Entry>& entries) {
        std::sort(entries.begin(), entries.end());

        SparseMatrix mat(r, c);
        if (entries.empty()) return mat;

        // Sum duplicates
        std::vector<Entry> compressed;
        if (!entries.empty()) {
            compressed.push_back(entries[0]);
            for (size_t i = 1; i < entries.size(); ++i) {
                if (entries[i].r == compressed.back().r && entries[i].c == compressed.back().c) {
                    compressed.back().v += entries[i].v;
                } else {
                    compressed.push_back(entries[i]);
                }
            }
        }

        mat.values.reserve(compressed.size());
        mat.col_indices.reserve(compressed.size());
        
        int current_row = 0;
        mat.row_ptr[0] = 0;
        
        for (const auto& e : compressed) {
            while (current_row < e.r) {
                mat.row_ptr[++current_row] = (int)mat.values.size();
            }
            mat.values.push_back(e.v);
            mat.col_indices.push_back(e.c);
        }
        
        while (current_row < r) {
            mat.row_ptr[++current_row] = (int)mat.values.size();
        }

        return mat;
    }

    /**
     * @brief Sparse Matrix-Vector product: y = A * x
     */
    std::vector<double> multiply(const std::vector<double>& x) const {
        if ((int)x.size() != cols) {
            throw std::runtime_error("SparseMatrix::multiply dimension mismatch.");
        }
        std::vector<double> y(rows, 0.0);
        for (int i = 0; i < rows; ++i) {
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                y[i] += values[k] * x[col_indices[k]];
            }
        }
        return y;
    }
};

} // namespace num
