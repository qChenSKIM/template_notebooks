##################################################
#  1. Tools > Terminal > New Terminal
#  Run in Unix Terminal Ctrl-Alt-Enter           #
##################################################
# we recommend running this is a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library("cmdstanr")
library(mvShapiroTest)
library(tidyverse)

# Set unix directories, /mnt/c/ at beginning for c:/
dir_data <- "1. CBC Recoding" # Where is R List with data
dir_model <- "99. Stan Code" # Location of Stan Code
dir_draws <- "2. HB Run/2. Stan HB Run/1. Output" # Where Stan stores draws.  Recommend a folder that does not sync
dir_out <- dir_draws # Location to put final output
# Load data and model
load(file.path(dir_data,"data_conjoint_RegularCoding.RData")) # Load data file
HB_model <- cmdstan_model(file.path(dir_model, "MNL_BartBlockCon_v1_6.stan"), quiet = TRUE, cpp_options = list(stan_threads = TRUE))
HB_model$print() # Just to verify

# Specify chains and threads
threads = list(parallel_chains = 2,
               threads_per_chain = 4)

# Specify constraints (sign only)
# For each parameter: 0 = no constraint, +1 = positive, -1 = negative
P <- data_conjoint$P
con_sign <- rep(0,P)
as.matrix(colnames(data_conjoint$code_master))
con_sign[21:39] <- -1 # Negative utilities for price slopes

# Modeling parameters. We include constraints above here.
# This overwrites/adds to the data file we pass to Stan.
data_model <- list(
  con_sign = con_sign,
  prior_cov = data_conjoint$prior_cov * 1, # Change cov scale here  
  df = 5, # Degrees of freedom
  prior_alpha = rep(0, P),
  a_sig = 10,
  cov_block = matrix(1, P, P),
  splitsize = round(.5 + data_conjoint$T/(4 * threads[[2]])),
  agg_model = NULL,
  tag = NULL
)

HB_model$sample(modifyList(data_conjoint, data_model),
                iter_warmup = 4,
                iter_sampling = 4,
                output_dir = dir_draws,
                chains = 2,
                parallel_chains = threads[[1]],
                threads_per_chain = threads[[2]],
                save_warmup = TRUE,
                refresh = 10,
                seed = 271,
                init = .1,
                show_messages = FALSE,
                validate_csv = FALSE
)

############################################################
##   Process Output
##########################################################
library("posterior")
csv_name <- c("MNL_BartBlockCon_v1_6-202107231212-1-498344.csv",
              "MNL_BartBlockCon_v1_6-202107231212-2-498344.csv"
) # You must specify names of output files in dir_draws

draws_upper <- read_cmdstan_csv(file.path(dir_draws, csv_name), variables = c("alpha"))
fit_stats <- summarize_draws(draws_upper$post_warmup_draws)

chain1_Alphas<- matrix(draws_upper$post_warmup_draws[,1,],
                       400,data_conjoint$P)

# Check normality
library('rstatix')
varChain<- as.data.frame(chain1_Alphas[,1])
shapTest<- shapiro_test(varChain[,1])

ggplot(varChain, aes(x=varChain[,1])) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.5, fill="#FF6666") +
      geom_vline(xintercept =c(quantile(varChain[,1],.05), quantile(varChain[,1],.95)), color = "red") + 
      geom_text(aes(x = quantile(varChain[,1],.005), y = 1, 
                    label = paste0("Shapiro p-value = ",round(shapTest$p.value,3))))


draws_beta <- read_cmdstan_csv(file.path(dir_draws, csv_name), variables = "beta_ind", sampler_diagnostics = "")
utilities <- matrix(colMeans(draws_beta$post_warmup_draws, dims = 2),
                    data_conjoint$I, data_conjoint$P, byrow = TRUE)
# Above assume draws were stored (P1,id1), (P2, id1), ..., (P1, id2)
# Otherwise byrow = FALSE
betas_final_r <- utilities %*% t(data_conjoint$code_master)
write.table(betas_final_r, file = file.path(dir_out, "betas_final_r.csv"), sep = ",", na = ".", row.names = FALSE)
write.table(fit_stats, file = file.path(dir_out, "fit_stats.csv"), sep = ",", na = ".", row.names = FALSE)

