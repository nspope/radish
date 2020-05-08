#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec backpropagate_laplacian_to_conductance (const arma::mat& tGl, const arma::mat& tGr, const arma::umat& tadj)
{
  if (tadj.n_rows != 2)
    Rcpp::stop("[backpropagate_gradient] tadj.n_rows != 2");

  if (tGl.n_rows != tGr.n_rows || tGl.n_cols != tGr.n_cols)
    Rcpp::stop("[backpropagate_gradient] dim(tGl) != dim(tGr)");

  const unsigned N = tGl.n_cols;

  arma::vec djj (N);
  for (unsigned vertex = 0; vertex < N; ++vertex)
    djj.at(vertex) = -arma::dot(tGl.col(vertex), tGr.col(vertex));

  arma::vec dl_dC (N+1, arma::fill::zeros);
  unsigned i, j;
  double dij;

  for (unsigned edge = 0; edge < tadj.n_cols; ++edge)
  {
    j = tadj.at(0, edge);
    i = tadj.at(1, edge);
    if (i <= j) 
      continue;
    else if (i > N)
      Rcpp::stop("[backpropagate_gradient] vertex out of bounds");
    else if (i == N)
    {
      dl_dC.at(i) += djj.at(j);
      dl_dC.at(j) += djj.at(j);
    }
    else
    {
      dij = djj.at(i) + djj.at(j) + 
        arma::dot(tGl.col(i), tGr.col(j)) +
        arma::dot(tGl.col(j), tGr.col(i));
      dl_dC.at(i) += dij;
      dl_dC.at(j) += dij;
    }
  }

  return dl_dC;
}

// [[Rcpp::export]]
arma::sp_mat backpropagate_conductance_to_laplacian (const arma::vec& dgrad__ddl_dC, const arma::umat& tadj)
{
  if (tadj.n_rows != 2)
    Rcpp::stop("[backpropagate_conductance] tadj.n_rows != 2");

  const unsigned N = dgrad__ddl_dC.n_elem - 1;

  unsigned i, j;
  double dij;

  arma::vec diagonal (N, arma::fill::zeros),
            offdiagonal (tadj.n_cols, arma::fill::zeros);

  for (unsigned edge = 0; edge < tadj.n_cols; ++edge)
  {
    j = tadj.at(0, edge);
    i = tadj.at(1, edge);

    if (i <= j)
      continue;
    else if (i > N)
      Rcpp::stop("[backpropagate_conductance] vertex out of bounds");

    dij = dgrad__ddl_dC.at(i) + dgrad__ddl_dC.at(j);
    diagonal.at(j) += dij;
    if (i != N)
    {
      offdiagonal.at(edge) = -dij;
      diagonal.at(i) += dij;
    }
  }

  arma::sp_mat dgrad__ddl_dQn (tadj, offdiagonal, N, N);
  dgrad__ddl_dQn.diag() = diagonal;

  return arma::trimatu(dgrad__ddl_dQn);
}
