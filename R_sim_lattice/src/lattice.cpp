
#include "misc_v12.h"
#include "probability_v16.h"

#include <Rcpp.h>

using namespace std;

//------------------------------------------------
// Simulate from lattice model assuming two alleles

// [[Rcpp::export]]
Rcpp::List sim_lattice_biallelic_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // get inputs from Rcpp format to base C++ format
  int demes_x = rcpp_to_int(args("demes_x"));
  int demes_y = rcpp_to_int(args("demes_y"));
  int N = rcpp_to_int(args("N"));
  double N_inv = 1.0 / double(N);
  double mu = rcpp_to_double(args("mu"));
  double m = rcpp_to_double(args("m"));
  vector<int> t_out = rcpp_to_vector_int(args("t_out"));
  int n_t_out = t_out.size();
  
  double p_init = rcpp_to_double(args("p_init"));
  bool torus = rcpp_to_bool(args["torus"]);
  bool silent = rcpp_to_bool(args["silent"]);
  
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  
  // objects for storing results
  Rcpp::List ret;
  
  // create two identical matrices for storing allele frequencies at all points
  // in a grid. We will step forward in time by alternating between these
  // matrices
  vector<vector<double>> p_mat1(demes_y, vector<double>(demes_x, p_init));
  vector<vector<double>> p_mat2(demes_y, vector<double>(demes_x));
  
  // option to store output at time 0
  int output_index = 0;
  int target_mat = 0;
  if (t_out[0] == 0) {
    ret.push_back(p_mat1);
    output_index++;
  }
  
  // loop through generations
  for (int t = 1; t < (max(t_out) + 1); ++t) {
    
    // report progress
    if (!silent) {
      update_progress(args_progress, "pb", t, max(t_out));
    }
    
    // update matrices
    if (target_mat == 0) {
      
      for (int i = 0; i < demes_y; ++i) {
        for (int j = 0; j < demes_x; ++j) {
          
          // deal with edges
          int up = i - 1;
          int down = i + 1;
          int left = j - 1;
          int right = j + 1;
          if (i == 0) {
            up = (torus) ? demes_y - 1 : i;
          }
          if (i == (demes_y - 1)) {
            down = (torus) ? 0 : i;
          }
          if (j == 0) {
            left = (torus) ? demes_x - 1 : j;
          }
          if (j == (demes_x - 1)) {
            right = (torus) ? 0 : j;
          }
          
          // update allele frequency
          double p_mig =  0.5*m*p_mat1[up][j] + 0.5*m*p_mat1[down][j] +
            0.5*m*p_mat1[i][left] + 0.5*m*p_mat1[i][right] + p_mat1[i][j]*(1 - 2*m);
          double p_mut = p_mig*(1 - mu*0.5) + (1 - p_mig)*mu*0.5;
          p_mat2[i][j] = rbinom1(N, p_mut) * N_inv;
        }
      }
      
    } else {
      
      for (int i = 0; i < demes_y; ++i) {
        for (int j = 0; j < demes_x; ++j) {
          
          // deal with edges
          int up = i - 1;
          int down = i + 1;
          int left = j - 1;
          int right = j + 1;
          if (i == 0) {
            up = (torus) ? demes_y - 1 : i;
          }
          if (i == (demes_y - 1)) {
            down = (torus) ? 0 : i;
          }
          if (j == 0) {
            left = (torus) ? demes_x - 1 : j;
          }
          if (j == (demes_x - 1)) {
            right = (torus) ? 0 : j;
          }
          
          // update allele frequency
          double p_mig =  0.5*m*p_mat2[up][j] + 0.5*m*p_mat2[down][j] +
            0.5*m*p_mat2[i][left] + 0.5*m*p_mat2[i][right] + p_mat2[i][j]*(1 - 2*m);
          double p_mut = p_mig*(1 - mu*0.5) + (1 - p_mig)*mu*0.5;
          p_mat1[i][j] = rbinom1(N, p_mut) * N_inv;
        }
      }
      
    }
    
    // save output
    if (t == t_out[output_index]) {
      if (target_mat == 0) {
        ret.push_back(p_mat1);
      } else {
        ret.push_back(p_mat2);
      }
      
      if (output_index < (n_t_out - 1)) {
        output_index++;
      }
    }
    
    // alternate target matrix
    target_mat = 1 - target_mat;
  }
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2 - t1);
  print("simulation completed in", time_span.count(), "seconds\n");
  
  // return Rcpp list
  return ret;
}

