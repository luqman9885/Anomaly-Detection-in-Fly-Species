library(MASS)
library(mvtnorm)
library(tidyverse)

setwd('../Dataset')

#Choose which cell compartment to use
#df <- read.csv('h10_dm_final_norm.csv')
df <- read.csv('h10_pa2r_final_norm.csv')

#Exclude albiceps mutant
df <- df[df$species != "Ch.albiceps_mutant",]
df$species <- factor(df$species)
df$family <- factor(df$family)
levels(df$family)

#Shuffle the rows so that it the species are randomized.
set.seed(17202806) 
df <- df[sample(nrow(df)), ]

#Exclude a1,b1,c1 from the first harmonic
df <- df[,-c(4,14,24)]

#Data Partitioning based on family levels
train_test_split <- function(FAMILY){
  test.data <- df[df$family == FAMILY,]
  test.data$label <- rep('unknown', nrow(test.data))
  train.data <- df[df$family != FAMILY,]
  
  test.data$family <- factor(test.data$family)
  train.data$family <- factor(train.data$family)
  
  mylist <- list()
  for(h in 1:5){
    test.set <- test.data
    train.set <- data.frame()
    
    for(i in 1:nlevels(train.data$family)){
      train.family <- train.data[train.data$family == levels(train.data$family)[i],]
      
      prop.20 <- floor(0.2*nrow(train.family))
      start <- 1 + prop.20*(h-1)
      if(h != 5){
        end <- prop.20 + prop.20*(h-1)
      }else{
        end <- nrow(train.family)
      }
      index <- start:end
      
      test.family <- train.family[index,]
      test.family$label <- rep('known', nrow(test.family))
      test.set <- rbind(test.set, test.family)
      
      train.family <- train.family[-index, ]
      train.family$label <- rep('known', nrow(train.family))
      train.set <- rbind(train.set, train.family)
    }
    mylist[[h]] <- list(train.set = train.set, test.set = test.set)
  }
  return(mylist)
}


result_list <- list()

confusion_family <- function(FAMILY){
  for(f in 1:5){
    train <- train_test_split(FAMILY)[[f]]$train.set
    test <- train_test_split(FAMILY)[[f]]$test.set
    
    #Train set
    train.info <- train[,c(1:3,41)]
    model <- lda(family ~. , data = train[,c(3,4:40)])
    train.lda <- predict(model)$x
    train.lda <- cbind(train.info,train.lda)
    
    #Estimated parameters for t-distribution
    d.f <- 3
    global.mean <- mean(train.lda$LD1)
    global.sd <- sd(train.lda$LD1)
    
    #Standardize LD1 for the t-distribution
    z <- (train.lda$LD1 - global.mean)/global.sd
    
    train.mat.result <- matrix(NA, nrow = nrow(train.lda), ncol = 2)
    
    for(j in 1:2){
      train.family <- train.lda[train.lda$family == levels(train.lda$family)[j],]
      mean.family <- mean(train.family$LD1)
      sd.family <- sd(train.family$LD1)
      
      log.known.density <- dnorm(train.lda$LD1, mean = mean.family, sd = sd.family,
                                 log = T)
      
      log.unknown.density <- dt(z, df = d.f, log = T) - log(global.sd)
      
      log.r <- log.known.density - log.unknown.density
      
      train.mat.result[,j] <- log.r
    }
    
    max.log.r <- apply(train.mat.result, 1, max)
    
    C <- quantile(max.log.r, 0.10)
    
    #Test set
    test.info <- test[,c(1:3,41)]
    test.lda <- predict(model, newdata = test[,c(3,4:40)])$x
    test.lda <- cbind(test.info,test.lda)
    
    test.mat.result <- matrix(NA, nrow(test.lda), ncol = 2)
    
    z_test <- (test.lda$LD1 - global.mean)/global.sd
    
    for(k in 1:2){
      train.family <- train.lda[train.lda$family == levels(train.lda$family)[k],]
      mean.family <- mean(train.family$LD1)
      sd.family <- sd(train.family$LD1)
      
      log.known.density <- dnorm(test.lda$LD1, mean = mean.family, sd = sd.family,
                                 log = T)
      
      log.unknown.density <- dt(z_test, df = d.f, log = T) - log(global.sd)
      
      log.r <- log.known.density - log.unknown.density
      
      test.mat.result[,k] <- log.r
    }
    
    max.log.r.test <- apply(test.mat.result,1,max)
    
    predicted.label <- numeric(nrow(test.lda))
    
    for(l in 1:nrow(test.lda)){
      if(max.log.r.test[l] <= C){
        predicted.label[l] <- 'unknown'
      }else{
        predicted.label[l] <- 'known'
      }
    }
    
    test.lda$predicted.label <- predicted.label
    predicted.index <- apply(test.mat.result,1,which.max)
    test.lda$predicted.family <- levels(train.lda$family)[predicted.index]
    
    result <- test.lda %>% 
      mutate(
        Actual.label = case_when(
          label == 'unknown' ~ 'unknown',
          TRUE ~ family
        ),
        Predict.label = case_when(
          predicted.label == 'unknown' ~ 'unknown',
          TRUE ~ predicted.family
        )
      ) %>% 
      select(family,species,Actual.label,Predict.label)
    
    family.labels <- unique(result$Actual.label)
    
    result_list[[f]] <- table(
      Predicted = factor(result$Predict.label, family.labels),
      Actual = factor(result$Actual.label, family.labels)
    )
  }
  return(result_list)
}

#Metrics Result average across five folds for one family

five_fold_result.family <- function(FAMILY){
  result <- confusion_family(FAMILY)
  
  five_fold_list <- list()
  
  for(f in 1:5){
    result1 <- result[[f]]
    
    sensitivity <- numeric(length = 3)
    
    for(i in 1:3){
      sensitivity[i] <- round(result1[i,i]/sum(result1[,i]),2)
    }
    
    names(sensitivity) <- colnames(result1)
    
    specificity <- numeric(length = 3)
    
    for(i in 1:3){
      specificity[i] <- round(sum(result1[-i,-i])/(sum(result1[i,-i]) + sum(result1[-i,-i])),2)
    }
    
    names(specificity) <- colnames(result1)
    
    proportion.unknown.predicted <- numeric(length = 3)
    for(i in 1:3){
      proportion.unknown.predicted[i] <- round(result1['unknown',i]/sum(result1[,i]),2)
    }
    
    names(proportion.unknown.predicted) <- colnames(result1)
    
    result.df <- data.frame(sensitivity,specificity,proportion.unknown.predicted)
    result.df$balanced.accuracy <- round((result.df$sensitivity + result.df$specificity)/2,2)
    result.df$class <- rownames(result.df)
    rownames(result.df) <- NULL
    
    five_fold_list[[f]] <- result.df
  }
  
  df.average <- data.frame()
  
  for(i in 1:5){
    df.average <- rbind(df.average,five_fold_list[[i]])
  }
  
  result <- df.average %>% 
    group_by(class) %>% 
    summarize(
      avg.sensitivity = round(mean(sensitivity),2),
      avg.specificity = round(mean(specificity),2),
      avg.balanced.accuracy = round(mean(balanced.accuracy),2),
      avg.proportion.unknown.predicted = round(mean(proportion.unknown.predicted),2)
    ) %>% 
    arrange(desc(avg.sensitivity))
  
  return(result)
}

#Results for 3 families
metrics_list <- list()

metrics_list <- lapply(1:3,function(x){
  five_fold_result.family(levels(df$family)[x])
  })

names(metrics_list) <- levels(df$family)

#Macro Result
macro.df <- data.frame()

for(i in 1:3){
  macro.df <- rbind(macro.df,metrics_list[[i]])
}

macro.result <- macro.df %>% 
  group_by(class) %>% 
  summarize(
    macro.sensitivity = round(mean(avg.sensitivity),2),
    macro.specificity = round(mean(avg.specificity),2),
    macro.balanced.accuracy = round(mean(avg.balanced.accuracy),2),
    macro.proportion.unknown.predicted = round(mean(avg.proportion.unknown.predicted),2)
  )

macro.result <- macro.result[-4,]
macro.result$capture.rate.unknown <- unlist(lapply(metrics_list, function(X) {X[X$class == 'unknown',2]}))

colnames(macro.result) <- c('family', 'sensitivity', 'specificity', 'balanced accuracy', 'false unknown rate', 'capture rate unknown')

#Summary performance of known families
macro.df <- macro.df[macro.df$class != 'unknown',]

performance.summary <- macro.df %>% 
  group_by(class) %>% 
  summarise(
    `macro sensitivity` = round(mean(avg.sensitivity),2),
    `se sensitivity` = round(sd(avg.sensitivity),2),
    `macro specificity` = round(mean(avg.specificity),2),
    `se specificity` = round(sd(avg.specificity),2),
    `macro balanced accuracy` = round(mean(avg.balanced.accuracy),2),
    `se balanced accuracy` = round(sd(avg.balanced.accuracy),2),
    `macro false unknown rate` = round(mean(avg.proportion.unknown.predicted),2),
    `se false unknown rate` = round(sd(avg.proportion.unknown.predicted),2)
  ) %>% 
  arrange(desc(`macro sensitivity`))

#Summary performance of unknown families
capture.list <- list()

capture.rate.family <- function(family){
  X <- confusion_family(family)
  
  capture.rate <- numeric(5)
  for(i in 1:5){
    conf <- X[[i]]
    capture.rate[i] <- conf['unknown','unknown']/sum(conf[,'unknown'])
  }
  
  mean.capture <- round(mean(capture.rate),2)
  se.capture <- round(sd(capture.rate),2)
  
  y <- c(mean.capture, se.capture)
  names(y) <- c('mean', 'SE')
  return(y)
} 

capture.result <- sapply(1:3, function(i){
  capture.rate.family(levels(df$family)[i])
})

colnames(capture.result) <- levels(df$family)

#Final Results 
#Performance of known families
performance.summary
#Performance of unknown families
t(capture.result)