######################
#      STAN TEST     #
#  PART 1: Set Stan  #
######################

################################
#            Step 1            #
#  Check if STAN is installed  #
################################

if(!require(cmdstanr)) {
  # Stan does not work well in version 3.6 of R
  if(grep("4.0",R.version$version.string) != 1){
    message("Error. Please download R version 4.0. Also make sure to have Rtools40.")
  }else{
    install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  }
  library(cmdstanr)
  install_cmdstan()
  cmdstan_path()
  
}

  #################
  # ADD LIBRARIES #
  #################
  
  if(!require(ellipsis)){
    install.packages("ellipsis")
  }
  library(ellipsis)
  
  if(!require(bayesplot)){
    install.packages("bayesplot")
  }
  library(bayesplot)
  
  if(!require(posterior)){
    install.packages("posterior")
  }
  library(posterior)
  
  if(!require(rsample)){
    install.packages("rsample")
  }
  library(rsample)


################################
#            Step 2            #
#        Set-up the data       #
################################

# Step 2 - Set-up the data
source("https://raw.githubusercontent.com/klattery/Estimation-Ecosystem/master/EE_Functions_Final.R")

### Specify Directories - with data and where to export results, Read Data into R ######
dir_data <- "C:/Users/n.camargo/SKIM/Swarovski - F4766 Swarovski Repositioning WTP/4. Analysis/1. CBC/CN/1. Necklace/1. CBC Recoding"
dir_out <- "C:/Users/n.camargo/SKIM/Swarovski - F4766 Swarovski Repositioning WTP/4. Analysis/1. CBC/CN/1. Necklace/1. CBC Recoding"

data_all <- read.csv(file.path(dir_data, "Recoded_CBC_CN_Necklace.csv"), as.is = TRUE)
`%notin%` <- Negate(`%in%`)

# Split data into training and test (test being a hold out task per respondent)

# Function to take holdout task
  #' @param data all data 
  #' @param idCol respondent ID col
  #' @param taskCol task col position
  #' @param numhold number of hold outs to use
select_hold <- function(data,idCol,taskCol, numhold = 1){
  ksplit <- split(data, data[,idCol])             # Split data by respondent
  set.seed(123)
  ksplit<- lapply(ksplit,function(x){
                         taskVec<- unique(x[,taskCol])                 # Check max tasks per respondent
                         holdOut<- sample(taskVec, numhold)            # Sample x random tasks
    
                         # Train Data
                         train_IDs<- x[,taskCol] %notin% holdOut       # Take train tasks
                         train_data<- x[train_IDs,]                    # Train data
                         # Test Data
                         test_IDs<- x[,taskCol] %in% holdOut           # Take Hold out tasks
                         test_data<- x[test_IDs,]                      # Test data
                         return(list(train_data = train_data, test_data = test_data))
  })
  
  train_data<- do.call(rbind,lapply(ksplit,'[[',"train_data"))
  test_data<- do.call(rbind,lapply(ksplit,'[[','test_data'))
  return(list(train_data = train_data, test_data = test_data))
}

splitData<- select_hold(data = data_all,
                        idCol = 1,
                        taskCol = 2,
                        numhold=1)

#train_data<- splitData$train_data
#test_data<- splitData$test_data

train_data<- data_all

#############################################
##          1.  Coding/ Data Prep          ##
#############################################
as.matrix(colnames(data_all))
idtaskdep <- train_data[,c(1,2,32)]  # Specify data with id, task, dep in that order  

indcode_spec <- list()
  # Code part-worths
  att_levels<- grep("Att",names(train_data))
    for(i in 1:length(att_levels)){
      indcode_spec[[i]]<- catcode(train_data,att_levels[i])
    }
  
  # Code user-code
  lengthList<- length(indcode_spec)
  slope_levels<- grep(c("Slope"),names(train_data))
  indcode_spec[[lengthList+1]]<- usercode(train_data,slope_levels)
  
  # Code none
  lengthList<- length(indcode_spec)
  none_levels<- grep(c("None"),names(train_data))
  indcode_spec[[lengthList+1]]<- usercode(train_data,none_levels)
  

indcode_list <- make_codefiles(indcode_spec)

# Check code_master is right
#View(indcode_list$code_master)
write.table(cbind(rownames(indcode_list$code_master), indcode_list$code_master),
            file = file.path(dir_out, "code_master.csv"),
            sep = ",", na = ".", row.names = FALSE)

#############################################
##    2.  Run Agg Model to Check Coding    ##
##        Save data files as R List        ##
#############################################
data_conjoint <- prep_file(idtaskdep, indcode_list)

# Add constraints
#constrain <- make_con(code_master=indcode_list$code_master,    # make_con(code_master) for no constraints
#                      con_specs = list(col_pos = NULL,
#                      col_neg = 30:54,
#                      row_rel = NULL),
#                      x0_try =  rep(0, data_conjoint$P)
#                      )
model_base <- list(
  func = list(
    pred = PredMNL,
    min = LL_Neg,
    gr = grad_MNL
  ),
  con = con_trivial(data_conjoint$P), # no constraints
  x0 =  rep(0, data_conjoint$P)       # start at 0
)

#model_base_con <- list(
#  func = list(
#    pred = PredMNL,
#    min = LL_Neg,
#    gr = grad_MNL
#  ),
#  con = as.matrix(constrain$constrain), # with constraints
#  x0 =  constrain$x0                    # start at 0
#)

agg_result <- agg_solve(data_conjoint, model_base) # Aggregate solution
betas <- data_conjoint$code_master %*% agg_result$par
#consistent<- FALSE                                     # For constraints

# For big data set this may take like 30 minutes
#JeffPriors <- JeffPrior(agg_result$par, data_conjoint, model_base,
#                        ObsFish = consistent, ExpFish = !consistent, score_n = 500)
#if (consistent){ 
#  cov_prior <- indprior %*% diag(1/JeffPriors$ObsFish)
#} else {
#  cov_prior <- data_conjoint$prior_cov %*% diag(1/JeffPriors$ExpFish)
#  data_conjoint$prior_cov<- cov_prior
#}


## Recommend reviewing betas to see if they look good
# If happy with agg model save coded data file
data_conjoint$tag <- "Regular Coding" # Specify a description/tag for your file
data_conjoint$agg_model$par_est <- agg_result$par
data_conjoint$agg_model$model_base <- model_base
save(data_conjoint, file = file.path(dir_out, "data_conjoint_RegularCoding.RData")) 

          