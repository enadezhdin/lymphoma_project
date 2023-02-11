###combined script for Lymphoma project
## author: E. Nadezhdin. en282@cam.ac.uk

library(tidyverse)
library(xtable)
library(RColorBrewer)
library(survival)
library(CoxHD)
library(Rcpp)
library(survAUC)
library(survminer)
library(randomForestSRC)
library(rpart)

## def functions
`%ni%` <- Negate(`%in%`)

training_idx_generator <- function(lymph_df, cols_exclude = c("PID", "OS_time", "OS_status", "PFS_time", "PFS_status")){
  total_idx <- c(1:dim(lymph_df)[2])
  idx_exclude <- which(colnames(lymph_df) %in% cols_exclude)
  train_set_idx <- total_idx[-idx_exclude]
  return(train_set_idx)
}

Lymph_OS_RFX <- function(coxRFX_os, data, x =365, tdMfBaseline = rep(1, ceiling(max(x))+1), tdMfAmlBaseline = rep(1, ceiling(max(x))+1), tdPrmBaseline = rep(1, ceiling(max(x))+1), tdOsBaseline = rep(1, ceiling(max(x))+1), ciType="analytical",  stage0="CP"){
  cppFunction('NumericVector computeHierarchicalSurvival(NumericVector x, NumericVector diffS0, NumericVector S1Static, NumericVector haz1TimeDep) {
                    int xLen = x.size();
                    double h;
                    NumericVector overallSurvival(xLen);
                    for(int i = 0; i < xLen; ++i) overallSurvival[i] = 1;
                    for(int j = 1; j < xLen; ++j) if(diffS0[j-1] != 0){
                    h = haz1TimeDep[j-1];
                    for(int i = j; i < xLen; ++i){
                    overallSurvival[i] += diffS0[j-1] * (1-pow(S1Static[i-j], h));
                    }
                    }
                    return overallSurvival;
                    }')
  ## Step 1: Compute KM survival curves and log hazard
  getS <- function(coxRFX, data, max.x=5000) {        
    if(!is.null(coxRFX$na.action)) coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
    data <- as.matrix(data[,match(colnames(coxRFX$Z),colnames(data)), drop=FALSE])
    r <- PredictRiskMissing(coxRFX, data, var="var2")
    H0 <- basehaz(coxRFX, centered = FALSE)
    hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
    x <- c(0:ceiling(max.x/25))*25/365
    S <- exp(-hazardDist(x))
    return(list(S=S, r=r, x=x, hazardDist=hazardDist, r0 = coxRFX$means %*% coef(coxRFX)))
  }
  
  kmOS <- getS(coxRFX = coxRFX_os, data = data, max.x=max(x))
  xx <- 0:ceiling(max(x)/25)
  sapply(1:nrow(data), function(i){
    os_Abs  <- cumsum(c(1,diff(kmOS$S^exp(kmOS$r[i,1]))))
    cbind(
      death_os = 1 - os_Abs,
      os_os = os_Abs,
      zeros = rep(0, 366) # this is dummy vector, to add third dimension to the second element of the array. Needed for plotting funct work correctly
    )
  }, simplify='array')
}

## plot function for individual predictions
plot_lymph_OS <-
  function(newdata,multistate, UID, i){
    lab=c(0,5,10,15,20,25) ## added
    pastel2<-c("#BFFFFF","#A3D8D1","#FFBFC9","#FF7387","#D1C299","#A3D8D1","#FF7387")
    newdata<-newdata[i,]; 
    sedimentPlot(-multistate[seq(1,361,30),1,i], x=seq(1,361,30),y0=1, y1=0,  col=c(pastel2[c(4,1)], "#D6E3DE"), xlab="", ylab="", xaxt="n",  main = as.character(UID))
    title(xlab="Time from diagnosis (years)",ylab="Proportion of patients", line = 2, cex.lab = 1) ###added
    lines(x=seq(1,361,30), y=1-rowSums(multistate[seq(1,361,30),c(1,3),i]), lwd=1) ### 
    segments(x0=0,y0=0.5,x1=newdata$OS_time/25,y1=0.5) ## defines line from beginning to death/or censor!
    points(x=newdata$OS_time/25,y=0.5,cex=newdata$OS_status*2,pch=19); ## defines black circle in case of death
    axis(side=1,at=lab*365.25/25,labels=lab) # added
    
  }

load(file = "sed_plot_function.RData")

set.seed(42)

## loading data

load("data/Lymphoma_full_dataset.RData")
lymph_work_ds <- lymph_full %>% dplyr::select(-cell_of_origin) %>% dplyr::filter(!is.na(IPI))

## plots on some discriptive stats
ages <- lymph_work_ds %>% group_by(age_gt_60) %>% summarize(n = n())
barplot(ages$n/length(lymph_work_ds$age_gt_60), main = "Patients age groups", ylim = c(0, 0.8), ylab = "% of patients", xlab = "age group", names.arg = c("age < 60", "age > 60"), cex.names = 1.5, cex.axis = 1.2, cex.lab = 1.2)

ipi_scores_groups <- lymph_work_ds %>% group_by(IPI) %>% summarize(n = n())
barplot(ipi_scores_groups$n/length(lymph_work_ds$IPI), main = "IPI scores distribution", ylab = "% of patients", xlab = "IPI score", ylim = c(0,0.4), names.arg = ipi_scores_groups$IPI, cex.names = 1.2, cex.lab = 1.2)

r_chop_surv <- Surv(lymph_work_ds$OS_time, lymph_work_ds$OS_status)
r_chop_fit <- survfit(r_chop_surv ~ rchop_treated, data = lymph_work_ds)
ggsurvplot(r_chop_fit, data = lymph_work_ds)

boxplot(lymph_work_ds$OS_time/365, main = "Distribution of Overall Survival \nfor the patients", ylab = "Time, years", xlab = "OS time", cex.names = 1.2, cex.axis = 1.2, cex.lab = 1.2)

genes_list <- colnames(lymph_work_ds)[(which(colnames(lymph_work_ds) == "PFS_status") + 1):length(colnames(lymph_work_ds))]
## creating df with only genet data

genet_lymph_df <- dplyr::select(lymph_work_ds, all_of(genes_list))

muts_per_ps <- rowSums(genet_lymph_df)

hist(muts_per_ps, breaks = 22, main = "Distribution of \nmutation counts per patient", xlim = c(0, 25), xlab = "# of mutations", ylab = "Counts")

mut_counts_df <- as.data.frame(apply(X = genet_lymph_df, 2, FUN = function(x) length(which(x == 1))))

mut_counts_df <- rownames_to_column(mut_counts_df)

colnames(mut_counts_df) <- c("genes", "counts")

mut_counts_df <- arrange(mut_counts_df, desc(counts))
mut_counts_df$genes <- str_split(mut_counts_df$genes, "_", simplify = T)[ ,1]
## plotting the df

ggplot(mut_counts_df) +
  geom_col(mapping = aes(x = reorder(genes, -counts), y = 100*(counts/679)), fill = 'blue') +  #ok, made monochrome
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 5, angle = 90)) +
  labs(x = "Genes", y = "% of patients")



##trainig a model on full dataset (w/o cell of origin column)

lymph_train_full_idx <- training_idx_generator(lymph_work_ds)

lymph_surv <- Surv(lymph_work_ds$OS_time, lymph_work_ds$OS_status)

lymph_full_model_fit <- CoxRFX(lymph_work_ds[ ,lymph_train_full_idx], lymph_surv, verbose = T)

lymph_model_preds <- Lymph_OS_RFX(lymph_full_model_fit, lymph_work_ds, x = 365*25)

survConcordance(lymph_surv ~ lymph_model_preds[20 * 365/25, 1, ])$concordance
# concordant 
# 0.7805156

### extracting weighted variances for the trained model haz coefficients
lymph_covariances <- cbind(colSums(cov(as.matrix(lymph_work_ds[ ,names(lymph_full_model_fit$coeff)]) %*% diag(lymph_full_model_fit$coefficients))), 
                  exp(lymph_full_model_fit$coefficients))

## formatting variances as a dataframe
model_variances <- tibble("Variable" = dimnames(lymph_covariances)[[1]], "Covariance" = lymph_covariances[ ,1], "LogHaz" = lymph_covariances[ ,2])
model_variances <- model_variances %>% arrange(desc(abs(Covariance)))

# adjusting genes names for plotting
model_variances$Variable[4] <- "TP53"
model_variances$Variable[5] <- "SGK1"
model_variances$Variable[6] <- "HIST1H1D"
model_variances$Variable[9] <- "CDKN2A"

barplot(model_variances$Covariance[1:10], main = "Model features covariances \ntop 10", ylim = c(0, 0.5), ylab = "Covariances", xlab = "Model feature", names.arg = model_variances$Variable[1:10], cex.names = 0.7, cex.axis = 1.2)

 

### performing k-fold cross validation of the model with split ratio = 0.8

split_ratio <- 0.8
total_idx <- c(1:dim(lymph_work_ds)[1])
concs_cross_val_test <- list("Train_conc" = 0, "Test_conc" = 0)

train_features_idx <- training_idx_generator(lymph_work_ds)


for (i in 1:50){
  running_train_set_idx <- sample(total_idx, as.integer(split_ratio*dim(lymph_work_ds)[1]))
  running_test_set_idx <- setdiff(total_idx, running_train_set_idx)
  running_surv <- Surv(lymph_work_ds$OS_time[running_train_set_idx], lymph_work_ds$OS_status[running_train_set_idx])
  running_model_fit <- CoxRFX(lymph_work_ds[running_train_set_idx, train_features_idx], running_surv)
  train_conc <- as.numeric(running_model_fit$concordance["concordance"])
  test_preds <- Lymph_OS_RFX(running_model_fit, lymph_work_ds[running_test_set_idx, ], x = 365*25)
  test_conc <- survConcordance(Surv(lymph_work_ds$OS_time[running_test_set_idx], lymph_work_ds$OS_status[running_test_set_idx]) ~ test_preds[5 * 365/25, 1, ])$concordance
  concs_cross_val_test$Train_conc[i] <- train_conc
  concs_cross_val_test$Test_conc[i] <- test_conc
}

## plotting the results
plot(concs_cross_val_test$Train_conc, ylim = c(0.6, 0.83), main = "Model k-fold cross-validation \n80/20 split", xlab = "training/testing iterations", ylab = "C-index", pch = 19, cex = 0.5, col = "red")
lines(concs_cross_val_test$Test_conc, col = "green")
legend(0, 0.67, legend = c("Training set", "Validation set"), lty=2:1, col = c("red", "green"), cex = 0.8)

boxplot(concs_cross_val_test$Train_conc, concs_cross_val_test$Test_conc, main = "Training/validation set C-indices \n50 iterations",
        names = c("Training set", "Validation set"), ylab = "C-index", cex.axis = 1.2)

mean(concs_cross_val_test$Train_conc) # 0.7845695

mean(concs_cross_val_test$Test_conc) # 0.7397962

## now running "-1 feature" model training

features_weights_test_fitting <- list("Feature" = 0, "Conc" = 0) #defing empty list to collect features and concordances

for (i in 0:length(train_features_idx)){
  if (i == 0){
    running_features_idx <- train_features_idx
    test_feature <- "full model"
  }else{
    running_features_idx <- train_features_idx[-i]
    test_feature <- colnames(lymph_work_ds)[train_features_idx[i]]
  }
  model_surv <- Surv(lymph_work_ds$OS_time, lymph_work_ds$OS_status)
  running_model_fit <- CoxRFX(lymph_work_ds[ , running_features_idx], model_surv)
  train_conc <- as.numeric(running_model_fit$concordance["concordance"])
  features_weights_test_fitting$Feature[i+1] <- test_feature
  features_weights_test_fitting$Conc[i+1] <- train_conc
}

## Plotting the results

plot(features_weights_test_fitting$Conc, pch = 19, cex = 0.3, col = "blue", ylab = "", xlab = "", ylim = c(0.75, 0.8), xaxt = "n")
title(main = "-1 feature mode training", xlab = "", ylab = "C-index")
axis(side = 1, at=c(1:121), labels = features_weights_test_fitting$Feature, las = 2, cex.axis = 0.5)

## formatting results as a dataframe.
model_concs_features <- as.data.frame(features_weights_test_fitting)
model_concs_features <- dplyr::arrange(model_concs_features, desc(Conc))

### training "minimal" model with only top 5 fetures based on the covariances values.

top_variance_features <- sort(abs(model_variances$Covariance), decreasing = T)[1:5]

top_var_model_idx <- which(colnames(lymph_work_ds) %in% names(top_variance_features))

### running 50 iters with min model

concs_cross_val_min_vars_test <- list("Train_conc" = 0, "Test_conc" = 0)

for (i in 1:50){
  running_train_set_idx <- sample(total_idx, as.integer(split_ratio*dim(lymph_work_ds)[1]))
  running_test_set_idx <- setdiff(total_idx, running_train_set_idx)
  running_surv <- Surv(lymph_work_ds$OS_time[running_train_set_idx], lymph_work_ds$OS_status[running_train_set_idx])
  running_model_fit <- CoxRFX(lymph_work_ds[running_train_set_idx, top_var_model_idx], running_surv)
  train_conc <- as.numeric(running_model_fit$concordance["concordance"])
  test_preds <- Lymph_OS_RFX(running_model_fit, lymph_work_ds[running_test_set_idx, ], x = 365*25)
  test_conc <- survConcordance(Surv(lymph_work_ds$OS_time[running_test_set_idx], lymph_work_ds$OS_status[running_test_set_idx]) ~ test_preds[5 * 365/25, 1, ])$concordance
  concs_cross_val_min_vars_test$Train_conc[i] <- train_conc
  concs_cross_val_min_vars_test$Test_conc[i] <- test_conc
}

plot(concs_cross_val_min_vars_test$Train_conc, ylim = c(0.6, 0.83), main = "RFX minimal model k-fold cross-validation", xlab = "training/testing iterations", ylab = "C-index", pch = 19, cex = 0.5, col = "red")
lines(concs_cross_val_min_vars_test$Test_conc, col = "blue")
legend(0, 0.67, legend = c("Training set", "Validation set"), lty=2:1, col = c("red", "blue"), cex = 0.8)

boxplot(concs_cross_val_min_vars_test$Train_conc, concs_cross_val_min_vars_test$Test_conc, main = "RFX minimal model training/validation set C-indices",
        names = c("Training set", "Validation set"), ylab = "C-index", cex.axis = 1.2)

mean(concs_cross_val_min_vars_test$Train_conc) # 0.7438802
mean(concs_cross_val_min_vars_test$Test_conc) # 0.7347456

### comparison of the RFX model performance with RandomForest non optimized model

## generating a new df to feed into rForest model

tmp_lymph_df <- lymph_work_ds[ ,c(2:10, 13:125)]
rForest_test <- rfsrc(Surv(OS_time, OS_status) ~ ., data = tmp_lymph_df, ntree = 100)

## calculating concordances for the rForest model
1 - get.cindex(time = tmp_lymph_df$OS_time, censoring = tmp_lymph_df$OS_status, predicted = rForest_test$predicted.oob) # 0.730429

## running k-fold cross-validations for rForest model
rforest_preds_cross_vals <- list("Train_conc" = 0, "Test_conc" = 0)
full_idx <- c(1:dim(tmp_lymph_df)[1])

for (i in 1:50){
  running_train_set_idx <- sample(full_idx, as.integer(split_ratio*dim(tmp_lymph_df)[1]))
  running_test_set_idx <- setdiff(full_idx, running_train_set_idx)
  running_model_fit <- rfsrc(Surv(OS_time, OS_status) ~ ., data = tmp_lymph_df[running_train_set_idx, ], ntree = 100)
  train_conc <- 1 - get.cindex(time = tmp_lymph_df$OS_time[running_train_set_idx], censoring = tmp_lymph_df$OS_status[running_train_set_idx], predicted = running_model_fit$predicted.oob)
  test_preds <- predict(running_model_fit, tmp_lymph_df[running_test_set_idx, c(1:7, 10:122)], importance = "none")
  test_conc <- 1 - get.cindex(time = tmp_lymph_df$OS_time[running_test_set_idx], censoring = tmp_lymph_df$OS_status[running_test_set_idx], predicted = test_preds$predicted)
  rforest_preds_cross_vals$Train_conc[i] <- train_conc
  rforest_preds_cross_vals$Test_conc[i] <- test_conc
}


plot(rforest_preds_cross_vals$Train_conc, ylim = c(0.6, 0.83), main = "rForest model k-fold cross-validation \n80/20 split", xlab = "training/testing iterations", ylab = "C-index", pch = 19, cex = 0.5, col = "blue")
lines(rforest_preds_cross_vals$Test_conc, col = "green")
legend(0, 0.67, legend = c("Training set", "Validation set"), lty=2:1, col = c("blue", "green"), cex = 0.8)

boxplot(rforest_preds_cross_vals$Train_conc, rforest_preds_cross_vals$Test_conc, main = "rForest model training/validation set C-indices \n50 iterations",
        names = c("Training set", "Validation set"), ylab = "C-index", cex.axis = 1.2)

## plotting res for two models
boxplot(concs_cross_val_test$Train_conc, concs_cross_val_test$Test_conc, rforest_preds_cross_vals$Train_conc, rforest_preds_cross_vals$Test_conc, main = "RFX and rForest models k-fold cross-validation",
        names = c("RFX model \ntraining", "RFX model \nvalidation", "rForest model \ntraining", "rForest model \nvalidation"), ylab = "C-index", cex.axis = 0.9)


### generatoing individual predictios plots based on the model preds outputs

for (i in 1:dim(lymph_work_ds)[1]){
  output_file_name <- sprintf("graphs_output/Lymph_test_%s.png", as.character(i))
  png(output_file_name)
  UID <- sprintf("patient_%s", as.character(lymph_work_ds$PID[i]))
  plot_lymph_OS(lymph_work_ds, lymph_model_preds, UID, i)
  dev.off()
}


### end of the script

