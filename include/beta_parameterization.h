/**
 * @file beta_parameterization.h
 * @brief C API for the Fortran beta parameterization library.
 *
 * Provides access to Legendre polynomial nuclear shape parameterization
 * from C, C++, Python (ctypes/cffi), and Rust.
 */

#ifndef BETA_PARAMETERIZATION_H
#define BETA_PARAMETERIZATION_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Compute nuclear radius grid from deformation parameters.
 *
 * Computes R(θ) on a uniform grid of n_grid points from θ=0 to θ=π
 * using the Legendre polynomial parameterization:
 *   R(θ) = R₀ [1 + Σ βλ Cλ Pλ(cosθ)]
 *
 * All arrays must be caller-allocated.
 *
 * @param[in]  params                Deformation parameters β₁..βₙ (array of n_params doubles)
 * @param[in]  n_params              Number of deformation parameters
 * @param[in]  n_grid                Number of θ grid points for output
 * @param[in]  n_quad                Number of Gauss-Legendre quadrature points (0 = use default)
 * @param[in]  apply_com_correction  Apply center-of-mass correction (0=false, nonzero=true)
 * @param[out] radii                 Output R(θ) values (array of n_grid doubles)
 * @param[out] corrected_beta10      The (possibly corrected) β₁₀ value
 * @param[out] is_valid              1 if shape is valid, 0 if invalid
 * @param[in]  message_buf_len       Size of message_buf including space for null terminator
 * @param[out] message_buf           Error message buffer (null-terminated on output)
 */
void beta_param_compute_radius_grid(
        const double *params, int n_params, int n_grid, int n_quad,
        int apply_com_correction,
        double *radii, double *corrected_beta10, int *is_valid,
        int message_buf_len, char *message_buf);

#ifdef __cplusplus
}
#endif

#endif /* BETA_PARAMETERIZATION_H */
