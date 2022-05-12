// ==============================================================================
// author          :Ghislain Vieilledent, Jeanne Clément
// email           :ghislain.vieilledent@cirad.fr, jeanne.clement16@laposte.net
// web             :https://ecology.ghislainv.fr
// license         :GPLv3
// ==============================================================================

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_erf.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

using namespace arma;
using namespace std;

//* ************************************************************ */
/* GSL mvgaussian */
int my_gsl_ran_multivariate_gaussian (const gsl_rng *r, const gsl_vector *mu, 
                                      const gsl_matrix *L, gsl_vector *result) {
  const size_t M = L->size1;
  const size_t N = L->size2;
  
  if (M != N) {
    GSL_ERROR("requires square matrix", GSL_ENOTSQR);
  }
  else if (mu->size != M) {
    GSL_ERROR("incompatible dimension of mean vector with variance-covariance matrix", GSL_EBADLEN);
  }
  else if (result->size != M) {
    GSL_ERROR("incompatible dimension of result vector", GSL_EBADLEN);
  }
  else {
    size_t i;
    for (i = 0; i < M; ++i) {
      gsl_vector_set(result, i, gsl_ran_ugaussian(r));
    }
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, L, result);
    gsl_vector_add(result, mu);
    return GSL_SUCCESS;
  }
}


/* ************************************************************ */
/* arma_mvgauss */

// [[Rcpp::export]]
arma::vec rmvnormRcpp(const gsl_rng *r, const arma::vec mu,
                      const arma::mat L) {
  
  // gsl vector mu
  gsl_vector *gsl_mu = gsl_vector_alloc(mu.n_elem);
  int size_mu = mu.n_elem;
  for (int i=0; i < size_mu; i++) {
    gsl_vector_set(gsl_mu, i, mu(i));
  }
  
  // gsl matrix L
  gsl_matrix *gsl_L = gsl_matrix_alloc(L.n_rows, L.n_cols);
  int nrows_L =  L.n_rows;
  int ncols_L =  L.n_cols;
  for (int i=0; i < nrows_L; i++) {
    for (int j=0; j < ncols_L; j++) {
      gsl_matrix_set(gsl_L, i, j, L(i,j));
    }
  }
  
  // gsl vector R
  gsl_vector *gsl_R = gsl_vector_alloc(mu.n_elem);
  gsl_vector_set_zero(gsl_R);
  
  // Call to gsl_ran_multivariate_gaussian
  my_gsl_ran_multivariate_gaussian(r, gsl_mu, gsl_L, gsl_R);
  
  // arma vec R
  arma::vec R; R.zeros(gsl_R->size);
  int size_R = gsl_R->size;
  for (int i=0; i < size_R; i++) {
    R(i) = gsl_vector_get(gsl_R, i);
  } 
  
  // free the memory
  gsl_vector_free(gsl_mu);
  gsl_matrix_free(gsl_L);
  gsl_vector_free(gsl_R);
  
  // Return result
  return R;
}

// End