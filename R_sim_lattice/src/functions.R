
#------------------------------------------------
# load other packages as needed
library(dplyr)
library(ggplot2)

#------------------------------------------------
# source more functions as needed for those included in this file
source("src/assertions_v14.R")
Rcpp::sourceCpp("src/lattice.cpp")

#------------------------------------------------
# update progress bar
# pb_list = list of progress bar objects
# name = name of this progress bar
# i = new value of bar
# max_i = max value of bar (close when reach this value)
# close = whether to close when reach end
#' @importFrom utils setTxtProgressBar
#' @noRd
update_progress <- function(pb_list, name, i, max_i, close = TRUE) {
  setTxtProgressBar(pb_list[[name]], i)
  if (i == max_i & close) {
    close(pb_list[[name]])
  }
}

#------------------------------------------------
#' @title Simulate allele frequencies from lattice model assuming two allele
#'   classes
#'
#' @description Simulate evolution under the lattice model of Kimura and Weiss.
#'   Can be defined over a plane or over a torus. Outputs the raw count of the
#'   reference allele in each deme, replicated over all specified output times
#'   as a list of matrices.
#'
#' @details Assumes a haploid population of size \code{N} and a single,
#'   biallelic locus. Implements a finite-alleles mutation model with equal
#'   chance of mutating to each allele. Initialises allele frequencies at
#'   \code{p_init} in all demes. Under the lattice model demes can only exchange
#'   migrants with other demes immediately adjacent, and do so with probability
#'   \eqn{m/2} in both the x- and y-dimensions. This means the total probability
#'   of migrating is actually \eqn{2m} for each individual in a deme, which
#'   seems strange but is how the lattice model was described in the original
#'   papers of Kimura and Weiss. By default, boundaries are reflecting, meaning
#'   demes at the edges have half the migration probability of demes further in.
#'   This tends to cause edge effects due to the buildup of relatedness at the
#'   boundaries. The option exists to make boundaries periodic, in which case
#'   the lattice model is defined over a torus rather than a plane. This
#'   eliminates edge effects, but leads to complex correlations over long
#'   distances.
#'   
#' @param demes_x,demes_y number of demes in each dimension.
#' @param N number of individuals per deme. The same for all demes.
#' @param mu mutation rate. Assumes finite-alleles model, with equal chance of
#'   mutating to each allele.
#' @param m per-generation probability of migrating to an adjacent deme in the
#'   x- and y-dimensions. The total probability of leaving a given deme is
#'   therefore \eqn{2m}.
#' @param t_out vector of times at which results will be output.
#' @param p_init the initial allele frequency, which is the same in all demes.
#' @param torus if \code{TRUE} then model is defined over a torus, i.e. with
#'   period rather than reflecting boundaries. Defaults to \code{FALSE}.
#' @param silent if \code{TRUE} then suppress output to console.
#'
#' @importFrom utils txtProgressBar
#' @export

sim_lattice_biallelic <- function(demes_x, demes_y, N, mu, m, t_out, p_init = 0.0,
                                  torus = FALSE, silent = FALSE) {
  
  # check inputs
  assert_single_pos_int(demes_x, zero_allowed = FALSE)
  assert_single_pos_int(demes_y, zero_allowed = FALSE)
  assert_single_pos_int(N, zero_allowed = FALSE)
  assert_single_bounded(mu)
  assert_single_bounded(m, right = 0.5)
  assert_vector_pos_int(t_out, zero_allowed = TRUE)
  assert_single_bounded(p_init)
  assert_single_logical(torus)
  assert_single_logical(silent)
  
  # make argument list
  args <- list(demes_x = demes_x, demes_y = demes_y, N = N, mu = mu, m = m,
               t_out = t_out, p_init = p_init, torus = torus, silent = silent)
  
  # create progress bars
  pb <- txtProgressBar(min = 0, max = max(c(1, t_out)), initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  # functions to pass to C++
  args_functions <- list(update_progress = update_progress)
  
  # run efficient C++ function
  output_raw <- sim_lattice_biallelic_cpp(args, args_functions, args_progress)
  
  #return(output_raw)
  
  # process output into matrices
  output_processed <- list()
  for (i in seq_along(t_out)) {
    output_processed[[i]] <- output_raw[[i]] %>%
      unlist() %>%
      matrix(nrow = demes_x)
  }
  names(output_processed) <- sprintf("t%s", seq_along(t_out))
  
  return(output_processed)
}

#------------------------------------------------
# plot matrix returned from sim_lattice_biallelic() as raster image using ggplot
plot_lattice_mat <- function(sim) {
  
  # get dimensions
  demes_x <- ncol(sim[[1]])
  demes_y <- nrow(sim[[1]])
  
  # get output into a dataframe, which is what ggplot expects
  df_freqs <- data.frame(sim = rep(seq_along(sim), each = demes_x*demes_y),
                         x = rep(1:demes_x, each = demes_y),
                         y = 1:demes_y,
                         p = unlist(sim))
  
  # plot raster
  df_freqs %>%
    ggplot() + theme_bw() +
    geom_raster(aes(x = x, y = y, fill = p)) +
    scale_fill_viridis_c(limits = c(0, 1), option = "magma") +
    facet_wrap(~sim)
}