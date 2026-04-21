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
    std::vector<double> values;     // CSR Values
    std::vector<int> col_indices;   // CSR Column Indices
    std::vector<int> row_ptr;       // CSR Row Pointers

    struct Entry {
        int r, c;
        double v;
        bool operator<(const Entry& other) const {
            if (r != other.r) return r < other.r;
            return c < other.c;
        }
    };

    std::vector<Entry> triplets; // Building buffer

    SparseMatrix(int r = 0, int c = 0) : rows(r), cols(c), row_ptr(r + 1, 0) {}

    void clear() {
        values.clear();
        col_indices.clear();
        row_ptr.assign(rows + 1, 0);
        triplets.clear();
    }

    /**
     * @brief Build CSR from the current triplets.
     */
    void compress() {
        if (triplets.empty()) return;
        std::sort(triplets.begin(), triplets.end());

        // Sum duplicates
        std::vector<Entry> compressed;
        compressed.push_back(triplets[0]);
        for (size_t i = 1; i < triplets.size(); ++i) {
            if (triplets[i].r == compressed.back().r && triplets[i].c == compressed.back().c) {
                compressed.back().v += triplets[i].v;
            } else {
                compressed.push_back(triplets[i]);
            }
        }

        values.clear();
        col_indices.clear();
        values.reserve(compressed.size());
        col_indices.reserve(compressed.size());
        
        row_ptr.assign(rows + 1, 0);
        int current_row = 0;
        row_ptr[0] = 0;
        
        for (const auto& e : compressed) {
            while (current_row < e.r) {
                row_ptr[++current_row] = (int)values.size();
            }
            values.push_back(e.v);
            col_indices.push_back(e.c);
        }
        
        while (current_row < rows) {
            row_ptr[++current_row] = (int)values.size();
        }
    }

    /**
     * @brief Convert to dense matrix (Vector of Vectors).
     */
    std::vector<std::vector<double>> to_dense() const {
        std::vector<std::vector<double>> dense(rows, std::vector<double>(cols, 0.0));
        for (int i = 0; i < rows; ++i) {
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                dense[i][col_indices[k]] = values[k];
            }
        }
        return dense;
    }

    /**
     * @brief Sparse Matrix-Vector product: y = A * x
     */
    std::vector<double> multiply(const std::vector<double>& x) const {
        if ((int)x.size() != cols) {
            throw std::runtime_error("SparseMatrix::multiply dimension mismatch.");
        }
        std::vector<double> y(rows, 0.0);
        #pragma omp parallel for
        for (int i = 0; i < rows; ++i) {
            double sum = 0;
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                sum += values[k] * x[col_indices[k]];
            }
            y[i] = sum;
        }
        return y;
    }
};

} // namespace num
