library(MASS)
library(mvtnorm)
library(tidyverse)

setwd('../Dataset')

#Choose which cell compartment to use 
df <- read.csv('h10_pa2r_final_norm.csv')
#df <- read.csv('h10_dm_final_norm.csv')

#Exclude albiceps mutant
df <- df[df$species != 'Ch.albiceps_mutant',]
df$species <- factor(df$species)
df$family <- factor(df$family)
levels(df$species)

#Exclude a1,b1,c1 from the first harmonic 
df <- df[,-c(4,14,24)]

#Data partitioning
train_test_split <- function(test_species){
  test.data <- df[df$species == test_species,]
  test.data$label <- rep('unknown',nrow(test.data))
  train.data <- df[df$species != test_species,]
  
  test.data$species <- factor(test.data$species)
  train.data$species <- factor(train.data$species)
  
  mylist <- list()
  for(h in 1:5){
    test.set <- test.data
    train.set <- data.frame()
    
    for(i in 1:nlevels(train.data$species)){
      train.sp <- train.data[train.data$species == levels(train.data$species)[i],]
      
      prop.20 <- floor(0.2*nrow(train.sp))
      start <- 1 + prop.20*(h-1)
      if(h != 5){
        end <- prop.20 + prop.20*(h-1)
      }else{
        end <- nrow(train.sp)
      }
      index <- start:end
      
      test.sp <- train.sp[index,]
      test.sp$label <- rep('known',nrow(test.sp))
      test.set <- rbind(test.set,test.sp)
      
      train.sp <- train.sp[-index,]
      train.sp$label <- rep('known',nrow(train.sp))
      train.set <- rbind(train.set,train.sp)
    }
    mylist[[h]] <- list(train.set = train.set, test.set = test.set)
  }
  return(mylist)
}

result_list <- list()

#Creating a function to perform the methodology on a particular species
result_species <- function(SPECIES, QUANTILE){
  
  for(f in 1:5){
    train <- train_test_split(SPECIES)[[f]]$train.set
    test <- train_test_split(SPECIES)[[f]]$test.set
    
    train.info <- train[,c(1:3,41)]
    
    model <- lda(species ~., data = train[,c(2,4:40)])
    train.lda <- predict(model)$x
    train.lda <- cbind(train.info,train.lda)
    
    d.f <- 3
    global_mean <- colMeans(train.lda[,5:15])
    global_cov <- cov(train.lda[,5:15])
    
    train.mat.result <- matrix(NA, nrow = nrow(train.lda), ncol = 12)
    
    #Calculating the threshold C using the training set 
    for(j in 1:12){
      train.species <- train.lda[train.lda$species == levels(train.lda$species)[j],]
      mean.species <- colMeans(train.species[,5:15])
      cov.species <- cov(train.species[,5:15])
      
      log.known.density <- dmvnorm(train.lda[,5:15], mean = mean.species,
                                   sigma = cov.species, log = T)
      
      log.unknown.density <- dmvt(train.lda[,5:15], delta = global_mean,
                                  sigma = global_cov, df = d.f, log = T)
      
      log.r <- log.known.density - log.unknown.density
      train.mat.result[,j] <- log.r
    }
    
    max.log.r <- apply(train.mat.result,1, max)
    
    C <- round(quantile(max.log.r,QUANTILE),2)
    
    #Testing on the test set 
    test.info <- test[,c(1:3,41)]
    test.lda <- predict(model, newdata = test[,c(2,4:40)])$x
    test.lda <- cbind(test.info,test.lda)
    
    test.mat.result <- matrix(NA, nrow = nrow(test.lda), ncol = 12)
    
    for(k in 1:12){
      train.species <- train.lda[train.lda$species == levels(train.lda$species)[k],]
      mean.species <- colMeans(train.species[,5:15])
      cov.species <- cov(train.species[,5:15])
      
      log.known.density <- dmvnorm(test.lda[,5:15], mean = mean.species,
                                   sigma = cov.species, log = T)
      
      log.unknown.density <- dmvt(test.lda[,5:15], delta = global_mean,
                                  sigma = global_cov, df = d.f, log = T)
      
      log.r <- log.known.density - log.unknown.density
      test.mat.result[,k] <- log.r
    }
    
    max.log.r.test <- apply(test.mat.result,1,max)
    
    predicted.label <- numeric(length = nrow(test.lda))
    
    for(l in 1:nrow(test.lda)){
      if(max.log.r.test[l] <= C){
        predicted.label[l] <- 'unknown'
      }else{
        predicted.label[l] <- 'known'
      }
    }
    
    test.lda$predicted.label <- predicted.label
    test.lda <- test.lda[,c(1,2,3,4,16,5:15)]
    
    predicted.index <- apply(test.mat.result,1,which.max)
    test.lda$predicted.species <- levels(train.lda$species)[predicted.index]
    
    test.lda <- test.lda[,c(1:5,17,6:16)]
    
    result <- test.lda %>% 
      mutate(
        Actual.label = case_when(
          label == 'unknown' ~ 'unknown',
          TRUE ~ species
        ),
        Predict.label = case_when(
          predicted.label == 'unknown' ~ 'unknown',
          TRUE ~ predicted.species
        )
      ) %>% 
      select(species,Actual.label,Predict.label)
    
    species.labels <- unique(result$Actual.label)
    
    result_list[[f]] <- table(
      Predicted = factor(result$Predict.label,species.labels), 
      Actual = factor(result$Actual.label,species.labels)
    )
  }
  
  return(result_list)
}

result.sp.quantile <- function(SPECIES,QUANTILE){
  result <- result_species(SPECIES,QUANTILE)
  
  five_fold_list <- list()
  
  for(f in 1:5){
    result1 <- result[[f]]
    
    sensitivity <- numeric(length = 13)
    
    for(i in 1:13){
      sensitivity[i] <- round(result1[i,i]/sum(result1[,i]),2)
    }
    
    names(sensitivity) <- colnames(result1)
    
    specificity <- numeric(length = 13)
    
    for(i in 1:13){
      specificity[i] <- round(sum(result1[-i,-i])/(sum(result1[i,-i]) + sum(result1[-i,-i])),2)
    }
    
    names(specificity) <- colnames(result1)
    
    proportion.unknown.predicted <- numeric(length = 13)
    for(i in 1:13){
      proportion.unknown.predicted[i] <- round(result1['unknown',i]/sum(result1[,i]),2)
    }
    
    names(proportion.unknown.predicted) <- colnames(result1)
    
    result.df <- data.frame(sensitivity,specificity,proportion.unknown.predicted)
    result.df$balanced.accuracy <- round((result.df$sensitivity + result.df$specificity)/2,2)
    result.df$class <- rownames(result.df)
    rownames(result.df) <- NULL
    result.df <- result.df[,c(5,1,2,4,3)]
    
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

metrics_list <- list()
metrics_list <- lapply(1:13, function(x){result.sp.quantile(levels(df$species)[x],0.10)})

names(metrics_list) <- levels(df$species)

macro.df <- data.frame()

for(i in 1:13){
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

macro.result <- macro.result[-14,]

macro.result$capture.rate.unknown <- unlist(lapply(metrics_list, function(X) {X[X$class == 'unknown',2]}))

family <- c(rep('Sarcophagidae',2),rep('Calliphoridae',7),
            rep('Sarcophagidae',2),'Muscidae','Sarcophagidae')

macro.result$family <- family
macro.result <- macro.result %>% arrange(desc(macro.sensitivity)) 
macro.df <- macro.df[macro.df$class != 'unknown',]

#Summary performance of known species 
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

#Summary performance of unknown species
capture.list <- list()

capture.rate.species <- function(species){
  X <- result_species(species,0.10)
  
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

capture.result <- sapply(1:13, function(i){
  capture.rate.species(levels(df$species)[i])
})

colnames(capture.result) <- levels(df$species)
#Final Results
#Performance of known species 
performance.summary
#Performance of unknown species
t(capture.result)