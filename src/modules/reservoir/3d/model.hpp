#pragma once
#include "lib/modules.hpp"
#include "state.hpp"
#include "lib/operators.hpp"
#include "modules/reservoir/well.hpp"
#include <cmath>
#include <vector>

namespace mod {
using namespace top;
namespace reservoir {

class Reservoir3DModel : public IModel {
private:
    double k;      // permeability [mD]
    double phi;    // porosity [fraction]
    double mu;     // viscosity [cP]
    double ct;     // total compressibility [psi^-1]
    double B;      // formation volume factor [rb/stb]
    
    std::vector<std::shared_ptr<ISourceSink>> sources;

    double diffusivity; 

public:
    Reservoir3DModel(double k_val, double phi_val, double mu_val, double ct_val, 
                     double B_val, const std::vector<std::shared_ptr<ISourceSink>>& sources_val)
        : k(k_val), phi(phi_val), mu(mu_val), ct(ct_val), B(B_val), 
          sources(sources_val) {
        
        diffusivity = 0.0002637 * k / (phi * mu * ct);
    }

    Vector evaluate_rhs(const IState& state) const override {
        const auto& r_state = dynamic_cast<const Reservoir3DState&>(state);
        int nx = r_state.spatial.nx, ny = r_state.spatial.ny, nz = r_state.spatial.nz;
        double dx = r_state.spatial.dx, dy = r_state.spatial.dy, dz = r_state.spatial.dz;
        
        Vector rhs(r_state.pressures.size(), 0.0);
        double inv_dx2 = 1.0 / (dx * dx);
        double inv_dy2 = 1.0 / (dy * dy);
        double inv_dz2 = 1.0 / (dz * dz);

        #pragma omp parallel for collapse(3)
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int cur = r_state.idx(i, j, k);
                    double p_cur = r_state.pressures[cur];
                    double p_w = (i == 0) ? p_cur : r_state.pressures[r_state.idx(i-1, j, k)];
                    double p_e = (i == nx - 1) ? p_cur : r_state.pressures[r_state.idx(i+1, j, k)];
                    double p_s = (j == 0) ? p_cur : r_state.pressures[r_state.idx(i, j-1, k)];
                    double p_n = (j == ny - 1) ? p_cur : r_state.pressures[r_state.idx(i, j+1, k)];
                    double p_b = (k == 0) ? p_cur : r_state.pressures[r_state.idx(i, j, k-1)];
                    double p_t = (k == nz - 1) ? p_cur : r_state.pressures[r_state.idx(i, j, k+1)];

                    double lap_p = (p_e - 2.0 * p_cur + p_w) * inv_dx2 + 
                                   (p_n - 2.0 * p_cur + p_s) * inv_dy2 + 
                                   (p_t - 2.0 * p_cur + p_b) * inv_dz2;
                    rhs[cur] = diffusivity * lap_p;
                }
            }
        }

        for (auto& s : sources) {
            s->apply(rhs, nullptr, r_state, 0.0, nullptr);
        }

        return rhs;
    }

    Vector build_residual(const IState& s_new, const IState& s_old, double dt) const override {
        const auto& r_new = dynamic_cast<const Reservoir3DState&>(s_new);
        const auto& r_old = dynamic_cast<const Reservoir3DState&>(s_old);
        int nx = r_new.spatial.nx, ny = r_new.spatial.ny, nz = r_new.spatial.nz;
        double dx = r_new.spatial.dx, dy = r_new.spatial.dy, dz = r_new.spatial.dz;
        
        Vector res(r_new.pressures.size());
        double tx = diffusivity * dt / (dx * dx);
        double ty = diffusivity * dt / (dy * dy);
        double tz = diffusivity * dt / (dz * dz);

        #pragma omp parallel for collapse(3)
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int cur = r_new.idx(i, j, k);
                    double p_cur = r_new.pressures[cur];
                    double p_w = (i == 0) ? p_cur : r_new.pressures[r_new.idx(i-1, j, k)];
                    double p_e = (i == nx - 1) ? p_cur : r_new.pressures[r_new.idx(i+1, j, k)];
                    double p_s = (j == 0) ? p_cur : r_new.pressures[r_new.idx(i, j-1, k)];
                    double p_n = (j == ny - 1) ? p_cur : r_new.pressures[r_new.idx(i, j+1, k)];
                    double p_b = (k == 0) ? p_cur : r_new.pressures[r_new.idx(i, j, k-1)];
                    double p_t = (k == nz - 1) ? p_cur : r_new.pressures[r_new.idx(i, j, k+1)];

                    double lap_p = (p_e - 2.0 * p_cur + p_w) * tx + 
                                   (p_n - 2.0 * p_cur + p_s) * ty + 
                                   (p_t - 2.0 * p_cur + p_b) * tz;
                    res[cur] = p_cur - r_old.pressures[cur] - lap_p;
                }
            }
        }

        Vector well_contributions(res.size(), 0.0);
        for (auto& s : sources) {
            s->apply(well_contributions, nullptr, r_new, dt);
        }
        #pragma omp parallel for
        for (int i = 0; i < (int)res.size(); ++i) {
            res[i] += dt * well_contributions[i]; 
        }

        return res;
    }

    SparseMatrix build_sparse_jacobian(const IState& state, double dt) const override {
        const auto& r_state = dynamic_cast<const Reservoir3DState&>(state);
        int nx = r_state.spatial.nx, ny = r_state.spatial.ny, nz = r_state.spatial.nz;
        double dx = r_state.spatial.dx, dy = r_state.spatial.dy, dz = r_state.spatial.dz;
        int n = (int)r_state.pressures.size();
        
        double tx = diffusivity * dt / (dx * dx);
        double ty = diffusivity * dt / (dy * dy);
        double tz = diffusivity * dt / (dz * dz);

        std::vector<SparseMatrix::Entry> entries;
        #pragma omp parallel
        {
            std::vector<SparseMatrix::Entry> local_entries;
            #pragma omp for
            for (int k = 0; k < nz; ++k) {
                for (int j = 0; j < ny; ++j) {
                    for (int i = 0; i < nx; ++i) {
                        int cur = r_state.idx(i, j, k);
                        local_entries.push_back({cur, cur, 1.0 + 2.0*tx + 2.0*ty + 2.0*tz});
                        
                        if (i > 0)      local_entries.push_back(SparseMatrix::Entry{cur, r_state.idx(i-1, j, k), -tx});
                        else            local_entries.back().v -= tx;
                        
                        if (i < nx - 1) local_entries.push_back(SparseMatrix::Entry{cur, r_state.idx(i+1, j, k), -tx});
                        else            local_entries.back().v -= tx;
                        
                        if (j > 0)      local_entries.push_back(SparseMatrix::Entry{cur, r_state.idx(i, j-1, k), -ty});
                        else            local_entries.back().v -= ty;
                        
                        if (j < ny - 1) local_entries.push_back(SparseMatrix::Entry{cur, r_state.idx(i, j+1, k), -ty});
                        else            local_entries.back().v -= ty;

                        if (k > 0)      local_entries.push_back(SparseMatrix::Entry{cur, r_state.idx(i, j, k-1), -tz});
                        else            local_entries.back().v -= tz;
                        
                        if (k < nz - 1) local_entries.push_back(SparseMatrix::Entry{cur, r_state.idx(i, j, k+1), -tz});
                        else            local_entries.back().v -= tz;
                    }
                }
            }
            #pragma omp critical
            {
                entries.insert(entries.end(), local_entries.begin(), local_entries.end());
            }
        }
        Vector dummy_res(n, 0.0);
        for (auto& s : sources) {
            s->apply(dummy_res, nullptr, r_state, dt, &entries);
        }
        return SparseMatrix::from_triplets(n, n, entries);
    }

    Matrix build_jacobian(const IState& state, double dt) const override {
        SparseMatrix S = build_sparse_jacobian(state, dt);
        Matrix J(S.rows, Vector(S.cols, 0.0));
        for (int i = 0; i < S.rows; ++i) {
            for (int k = S.row_ptr[i]; k < S.row_ptr[i + 1]; ++k) {
                J[i][S.col_indices[k]] = S.values[k];
            }
        }
        return J;
    }

    Vector apply_jacobian(const IState& state, const Vector& v, double dt) const override {
        const auto& r_state = dynamic_cast<const Reservoir3DState&>(state);
        int nx = r_state.spatial.nx, ny = r_state.spatial.ny, nz = r_state.spatial.nz;
        double dx = r_state.spatial.dx, dy = r_state.spatial.dy, dz = r_state.spatial.dz;
        
        double tx = diffusivity * dt / (dx * dx);
        double ty = diffusivity * dt / (dy * dy);
        double tz = diffusivity * dt / (dz * dz);

        Vector res(v.size());
        #pragma omp parallel for collapse(3)
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    int cur = r_state.idx(i, j, k);
                    
                    // (I - dt * eta * L) * v
                    double v_cur = v[cur];
                    double v_w = (i == 0) ? v_cur : v[r_state.idx(i-1, j, k)];
                    double v_e = (i == nx - 1) ? v_cur : v[r_state.idx(i+1, j, k)];
                    double v_s = (j == 0) ? v_cur : v[r_state.idx(i, j-1, k)];
                    double v_n = (j == ny - 1) ? v_cur : v[r_state.idx(i, j+1, k)];
                    double v_b = (k == 0) ? v_cur : v[r_state.idx(i, j, k-1)];
                    double v_t = (k == nz - 1) ? v_cur : v[r_state.idx(i, j, k+1)];
                    
                    double d2v_dx2_dt_eta = (v_e - 2.0 * v_cur + v_w) * tx;
                    double d2v_dy2_dt_eta = (v_n - 2.0 * v_cur + v_s) * ty;
                    double d2v_dz2_dt_eta = (v_t - 2.0 * v_cur + v_b) * tz;
                    
                    res[cur] = v_cur - (d2v_dx2_dt_eta + d2v_dy2_dt_eta + d2v_dz2_dt_eta);
                }
            }
        }
        return res;
    }
};

} // namespace reservoir
} // namespace mod
