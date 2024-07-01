# example_lattice.R
#
# Author: Bob Verity
# Date: 2023-03-10
#
# Purpose:
# Simple deploy script showing how to simulate from the lattice model and plot
# results.
#
# ------------------------------------------------------------------

# source R functions
source("src/functions.R")

# define parameters
demes_x <- 32
demes_y <- 32
N <- 1e3
mu <- 1e-5
m <- 5*1e-3
t_out <- seq(0, 5e4, l = 6)


# simulate from lattice model
sim1 <- sim_lattice_biallelic(demes_x = demes_x, demes_y = demes_y,
                              N = N, mu = mu, m = m, t_out = t_out,
                              torus = TRUE)

# get average allele frequency over the entire matrix at each point in time, and
# plot. This can be used to work out how long it takes for the process to reach
# equilibrium (typically about 5/mu generations)
df_mean <- data.frame(t = t_out,
                      p_mean = mapply(mean, sim1))
df_mean %>%
  ggplot() + theme_bw() +
  geom_line(aes(x = t, y = p_mean)) +
  ylim(c(0, 1)) +
  xlab("Time") + ylab("Mean allele frequency")

# plot matrices to visualise spatial autocorrelation
plot_lattice_mat(sim1)


n_sims <- 10000

for (i in 1:n_sims){
  print(i)
  
  N_df  <- rep(N, length(t_out))
  mu_df <- rep(mu, length(t_out))
  m_df  <- rep(m, length(t_out))
  
  df_params <- data.frame(cbind(N_df, mu_df, m_df))
  names(df_params) <- c('N', 'mu', 'm')
  
  # simulate from lattice model
  sim1 <- sim_lattice_biallelic(demes_x = demes_x, demes_y = demes_y,
                                N = N, mu = mu, m = m, t_out = t_out,
                                torus = TRUE)
  
  df_sim <- data.frame(matrix(unlist(sim1), nrow=length(sim1), byrow=TRUE))
  names(df_sim) <- paste0('p', c(1: (demes_x * demes_y)))
  
  df_sim <- cbind(df_params, df_sim)
  
  if (i==1){
    df <- df_sim
  } else {
    df <- rbind(df, df_sim)
  }
  if(i %% 1000 == 0) {
	fname <- paste0("/home/scratch/popgen/samples/samps_N",  N, "_mu", mu, "_m", m, ".csv")
	write.csv(df, fname, row.names=FALSE)
  }
}



