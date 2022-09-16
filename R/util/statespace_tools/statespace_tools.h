#ifndef STATESPACE_TOOLS_H
#define STATESPACE_TOOLS_H

/**
 * export for use via nimble::nimbleExternalCall in R
 *
 * nimble requires all non-scalar model data be passed as double arrays, and
 * the supporting functions accommodate nimble's requirement.
 */
double nimLLCompressedRaw(double* obs_lik_dict, double* obs, double* txmat_dict,
                          double* txmat_seq, double* x0, int nstates, int nt);

/**
 * export for use via nimble::nimbleExternalCall in R
 *
 * nimble requires all non-scalar model data be passed as double arrays, and
 * the supporting functions accommodate nimble's requirement.
 */
double nimLL2LayerCompressedRaw(double* obs_lik_dict, double* obs,
                                double* txmat_dict, double* txmat_seq,
                                double* x0, int num_obs_states,
                                int num_latent_states, int nt);

/**
 * export for use via nimble::nimbleExternalCall in R
 *
 * nimble requires all non-scalar model data be passed as double arrays, and
 * the supporting functions accommodate nimble's requirement.
 */
double nimLL2LayerPartialRaw(double* obs_lik_dict, double* obs,
                             double* txmat_seq, double* x0, int num_obs_states,
                             int num_latent_states, int nt, bool logEntries);

#endif
