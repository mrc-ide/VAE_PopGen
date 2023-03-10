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
mu <- 1e-4
m <- 1e-2
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
