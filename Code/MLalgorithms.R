#### Dravet Prediction Model for predicting Severity of Epilepsy

##Objective: Use Claims Level Information to predict Severity of Epilepsy

#Overview for Predictive Modeling####
#0. Exploration
#1. Data Pre-Processing
#2. Feature Selection
#2. Data splitting (train/test)
#3. Model training and tuning
#4. Performance evaluation

#0. Introduction####


#Objectives
#1. Load libraries and the data
#2. Transform data for computational efficiency
#3. Remove variables that are unnecessary (i.e. highly corelated or not predictive)
rm(list=ls())
#Load libraries and data
install.packages("corrplot")
install.packages("AppliedPredictiveModeling")

library(caret)
library(ggplot2)
library(plyr)
library(reshape2)
library(corrplot)
library(AppliedPredictiveModeling)
transparentTheme(trans = .4)

#Comorbidity data
comorb <- read.csv("C:/Users/ck/Dropbox/Research/Jon and CK/Dravet/Aim3 Algorithm/Data/Complex Chronic Conditions.csv")
comorb <- comorb[,c(1:2,397:410)] #Only the last 13 columns that summarize the chronic conditions
#Dravet CLaims data
drav <- read.csv("C:/Users/ck/Dropbox/Research/Jon and CK/Dravet/Aim3 Algorithm/Data/Cases and Controls_Summary Dataset_One Record per Patient_11-21-2016.csv")
#Medication data
med <- read.csv("C:/Users/ck/Dropbox/Research/Jon and CK/Dravet/Aim3 Algorithm/Data/medications.csv")

#Use top 5 dravet medications
med$clobazam <- ifelse(med$GENERIC_NAME == "CloBAZam", 1, 0)
med$midazolam <- ifelse(med$GENERIC_NAME == "Midazolam HCl", 1, 0)
med$DIAZEPAM <- ifelse(med$GENERIC_NAME == "DiazePAM", 1, 0)
med$Divalproex <- ifelse(med$GENERIC_NAME == "Divalproex Sodium", 1, 0)
med$LevETIRAcetam <- ifelse(med$GENERIC_NAME == "LevETIRAcetam", 1, 0)

#Top 1 for CAE
med$CLONAZEPAM <- ifelse(med$GENERIC_NAME == "ClonazePAM", 1, 0)

#Top 1 for ISLG
med$nacl <- ifelse(med$GENERIC_NAME == "Sodium Chloride", 1, 0)

#Make one patient one record for medication

t <- med[which(med$GENERIC_NAME=="CloBAZam" | med$GENERIC_NAME == "Midazolam HCl" | med$GENERIC_NAME == "DiazePAM" | med$GENERIC_NAME == "Divalproex Sodium" | med$GENERIC_NAME == "LevETIRAcetam" | med$GENERIC_NAME == "ClonazePAM" | med$GENERIC_NAME == "Sodium Chloride"),]
l<-count(t, c("UNIQUE_ID","GENERIC_NAME"))

w <- reshape(l, 
             timevar = "GENERIC_NAME",
             idvar = c("UNIQUE_ID"),
             direction = "wide")
w <- as.data.frame(w)
w[is.na(w)] <- 0

total <- merge(comorb,drav,by="UNIQUE_ID")

write.csv(w, file = "C:/Users/ck/Dropbox/Research/Jon and CK/Dravet/Aim3 Algorithm/Data/med_freq.csv")
med_freq <- read.csv("C:/Users/ck/Dropbox/Research/Jon and CK/Dravet/Aim3 Algorithm/Data/med_freq.csv", row.names = 1)

#Remove unnecessary data file
l <- NULL
med <- NULL
t <- NULL
w <- NULL

temp <- merge(total, med_freq, by="UNIQUE_ID")

#Visualize age, total_charges, Prescription_count, and freq.DiazePAM
myvars <- names(final) %in% c("AGE_YEARS", "TOTAL_CHARGES","PRESCRIPTION_COUNT", "freq.DiazePAM") 
featurePlot(x = final[,myvars], 
            y = as.factor(final$Class), 
            plot = "pairs",
            ## Add a key at the top
            auto.key = list(columns = 2))

#Overlayed density plot
transparentTheme(trans = .9)
featurePlot(x = final[myvars], 
            y = as.factor(final$Class),
            plot = "density", 
            ## Pass in options to xyplot() to 
            ## make it prettier
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")), 
            adjust = 1.5, 
            pch = "|", 
            layout = c(4, 1), 
            auto.key = list(columns = 3))
myvars <- names(final) %in% c("freq.CloBAZam","freq.Divalproex.Sodium", "freq.LevETIRAcetam", "freq.ClonazePAM") 
featurePlot(x = final[myvars], 
            y = as.factor(final$Class),
            plot = "density", 
            ## Pass in options to xyplot() to 
            ## make it prettier
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")), 
            adjust = 1.5, 
            pch = "|", 
            layout = c(4, 1), 
            auto.key = list(columns = 3))

#1. Data management####


#Preprocess data

#Make vectors with id based on the dataset for future use
comorb_ID <- comorb$UNIQUE_ID
drav_ID <- drav$UNIQUE_ID
med_ID <- med_freq$UNIQUE_ID

#Make data with only conntinuous variables for center,scale,and potentially PCA
comorb_proc <- comorb[,c(1,15)] #comorb data with only continuous variable
myvars <- names(drav) %in% c("COHORT_ID", "MATCHED_TRIAD","INDEX_DATE_SHIFT","DISTANCE","MED_GIVEN_OR_PRESCRIBED","X") 
drav_proc <- drav[!myvars] #drav data with only continuous variables

pre_proc <- merge(comorb_proc,drav_proc,by="UNIQUE_ID")
pre_proc <- merge(pre_proc, med_freq, by="UNIQUE_ID")

#Remove Unique ID for later merging
pre_proc_ID <- pre_proc$UNIQUE_ID 
pre_proc_SEX_NUM <- pre_proc$SEX_NUM
pre_proc_INSURANCE <- pre_proc$INSURANCE
myvars <- names(pre_proc) %in% c("UNIQUE_ID", "SEX_NUM","INSURANCE") 
pre_proc <- pre_proc[!myvars]

#BoxCox transform, center, scale, and PCA

trans <- preProcess(pre_proc, method=c("BoxCox","center","scale","nzv"))
trans
transformed <- predict(trans,pre_proc)


#merge back the unique id for merging the ccc to make dummy variables
transformed$UNIQUE_ID <- pre_proc_ID
transformed$SEX_NUM <- pre_proc_SEX_NUM
transformed$INSURANCE <- pre_proc_INSURANCE
transformed <- merge(transformed, comorb[,1:14], by="UNIQUE_ID")


#Change variables that are factor variables
names <- c('neuromusc_ccc','cvd_ccc','respiratory_ccc','congeni_genetic_ccc','tech_dep_ccc',"renal_ccc","GI_ccc","hemato_immu_ccc","metabolic_ccc","malignancy_ccc","neonatal_ccc","transplant_ccc",'SEX_NUM', 'INSURANCE')
transformed[,names] <- lapply(transformed[,names] , factor)

#Make integers into numeric now and change COHORT to integer afterwards as outcome
transformed[sapply(transformed, is.integer)] <- lapply(transformed[sapply(transformed, is.integer)], as.numeric)
transformed$COHORT <- as.integer(transformed$COHORT)

#Dummyfy the factors
transformedDummy <- dummyVars("~.",data=transformed, fullRank=F)
transformed <- as.data.frame(predict(transformedDummy,transformed))

#Rename COHORT to outcome
outcomeName <- 'COHORT'
predictorsNames <- names(transformed)[names(transformed) != outcomeName]

#Missing data check
transformed <- na.omit(transformed)

#Combine groups into severe vs mild
transformed$COHORT2 <- ifelse(transformed$COHORT==1,'Nonsevere','Severe')
transformed$COHORT2 <- as.factor(transformed$COHORT2)
outcomeName <- 'COHORT2'

#Taking out nearzerovariance variables finally
nzv <- nearZeroVar(transformed)
filteredtotal <- transformed[, -nzv]

#Optional step to check correlations although we removed everything by using PCA
#Check correlation of continuous variables
correlations <- cor(pre_proc)
corrplot(correlations, order="hclust")

#Find correlations that are high 0.75
highCorr <- findCorrelation(correlations, cutoff=.75)
length(highCorr)

#Make final before feature selection
myvars <- names(filteredtotal) %in% c("UNIQUE_ID", "COHORT") 
final <- filteredtotal[,!myvars]
final$Class <- final$COHORT2
final$COHORT2 <- NULL


#3. Data Splitting####

#Objective: Undergo test and train data splitting for determining optimal tuning parameter and model

#Set seed to reproduce the data
set.seed(1)

trainingrow <- createDataPartition(
  filteredtotal$COHORT2,
  p=0.8,
  list=FALSE
)

head(trainingrow)

#Subset data into predictors and classes
myvars <- names(filteredtotal) %in% c("UNIQUE_ID", "COHORT2","COHORT") 
trainPredictors <-  filteredtotal[trainingrow, !myvars]
myvars <- names(filteredtotal) %in% c( "COHORT2") 
trainClasses <- filteredtotal[trainingrow, myvars]

#create test set
myvars <- names(filteredtotal) %in% c("UNIQUE_ID", "COHORT2","COHORT") 
testPredictors <-  filteredtotal[-trainingrow, !myvars]
myvars <- names(filteredtotal) %in% c( "COHORT2") 
testClasses <- filteredtotal[-trainingrow, myvars]


#Resampling test/train 3 times

repeatedSplits <- createDataPartition(
  trainClasses,
  p=0.8,
  times = 3
)

# 10 fold cross validation
cvSplits <- createFolds(trainClasses,
                        k=10,
                        returnTrain = TRUE)

str(cvSplits)
fold1 <- cvSplits[[1]] # get the first fold
cvPredictors1 <- trainPredictors[fold1,]
cvClasses1 <- trainClasses[fold1]

nrow(trainPredictors)
nrow(cvPredictors1) # There is a difference bc cross validation resampling doesn't use the whole sample


#Use the final dataset that contains pre-processed predictors and the cohort information
set.seed(1234)
inTrain <- createDataPartition(final$Class, p = .65)[[1]]
finalTrain <- final[ inTrain, ]
finalTest  <- final[-inTrain, ]

#Check proportion of severe vs. nonsevere to see if data is imbalanced
prop.table(table(finalTest$Class)) #nonsevere only 22% so oversample nonsevere

#install.packages("DMwR")
library(DMwR)

finalTrain <- SMOTE(Class ~ ., finalTrain, k=5 ,perc.over = 500, perc.under=100)
prop.table(table(finalTrain$Class)) #nonsevere only 22% so oversample nonsevere
prop.table(table(finalTest$Class)) #nonsevere only 22% so oversample nonsevere

#4. Model buidling####

#We will be fitting 4 models
#CART, Boosted Tree, SVM, and GLM

#4.1 CART MODEL#####

library(rpart)
#rpart1 <- rpart(Class ~ ., data = finalTrain,
#                control = rpart.control(maxdepth = 5))
#rpart1

#recursive partitioning and plotting
#install.packages("partykit")
library(partykit)
#rpart1a <- as.party(rpart1)

#plot(rpart1a)

rpartfull <- rpart(Class ~ . , data=finalTrain)
rpartfull
plot(rpartfull)
partPred<-predict(rpartfull,finalTest,type="class")
confusionMatrix(partPred, finalTest$Class)

#Manual tuning using caret package

cvCtrl <- trainControl(method = "repeatedcv", repeats = 3,
                       summaryFunction = twoClassSummary
                       ,
                       classProbs = TRUE,
                       sampling = "up"
)
set.seed(123)
rpartcar <- train(Class ~ . , data=finalTrain, method = "rpart",
                  tuneLength= 30,
                  metric = "Sens",
                  trControl = cvCtrl )

rpartcar
plot(rpartcar)

#Predict classess using rpart
rpartPred2<-predict(rpartcar,finalTest)
confusionMatrix(rpartPred2,finalTest$Class)
rpartProbs<-predict(rpartcar,finalTest,type="prob")

rpartROC<-roc(finalTest$Class,rpartProbs[,"Severe"])
plot(rpartROC,type="S",print.thres=.5, xlim=c(1,0))
rpartROC


#4.2 Random forest#####
library(randomForest)
set.seed(123)
rf <- train(Class ~ . , data=finalTrain, method = "rf",
            metric = 'Sens',
            trControl = cvCtrl )
plot(rf)
rf #High specificity but low sensitivity

#Predict classess using rf
rfPred <-predict(rf,finalTest)
confusionMatrix(rfPred,finalTest$Class)
rfProbs<-predict(rf,finalTest,type="prob")
rfROC <- roc(finalTest$Class,rfProbs[,"Severe"])



#AUC values

rfROC
rpartROC


#4.3 Boosted Tree algorithm don't use#####
#install.packages("C50")
library(C50)

grid <- expand.grid(.model = "tree",
                    .trials = c(1:100),
                    .winnow = FALSE)

c5Tune <- train(finalTrain, finalTrain$Class,
                method = "C5.0"
                ,
                metric = "ROC",
                tuneGrid = grid
                ,
                trControl = cvCtrl)
c5Tune
plot(c5Tune)

c5Pred <-predict(c5Tune,finalTest)
confusionMatrix(c5Pred,finalTest$Class)

c5Probs<-predict(c5Tune,finalTest,type="prob")

c5ROC <- roc(predictor =  c5Probs$Severe,
             response = finalTest$Class)
c5ROC
plot(c5ROC,type="S")

#4.4 SVM#####

## The model fitting code shown in the computing section is fairly
## simplistic.  For the text we estimate the tuning parameter grid
## up-front and pass it in explicitly. This generally is not needed,
## but was used here so that we could trim the cost values to a
## presentable range and to re-use later with different resampling
## methods.
#install.packages("kernlab")
library(kernlab)
set.seed(231)
sigDist <- sigest(Class ~ ., data = finalTrain, frac = 1)
svmTuneGrid <- expand.grid(sigma = sigDist, C = 2^(-2:7))

### Optional: parallel processing can be used via the 'do' packages,
### such as doMC, doMPI etc. We used doMC (not on Windows) to speed
### up the computations.

### WARNING: Be aware of how much memory is needed to parallel
### process. It can very quickly overwhelm the available hardware. We
### estimate the memory usage (VSIZE = total memory size) to be 
### 2566M/core.

set.seed(1056)
svmFit <- train(Class ~ .,
                data = final,
                method = "svmRadial",
                tuneGrid = svmTuneGrid,
                trControl = trainControl(method = "repeatedcv", 
                                         repeats = 5,
                                         classProbs = TRUE,
                                         selectionFunction = "tolerance"))

set.seed(1056)
svmFit <- train(Class ~ .,
                data = finalTrain,
                method = "svmRadial",
                metric = 'Sens',
                tuneGrid = svmTuneGrid,
                trControl = cvCtrl)


library(doParallel) #SVM with recursive feature extraction... takes too long
registerDoParallel(cores=3)
svmFit <- rfe(Class ~ .,
                data = finalTrain,
                sizes = c(2,5,10,30),
                rfeControl = rfeControl(functions = caretFuncs,
                                        number = 50), #number indicates recursive outer sampling number
                method = "svmRadial",
                tuneGrid = svmTuneGrid,
                trControl = cvCtrl)

svmFit$finalModel

plot(svmFit, metric = "ROC",
     scales = list(x = list(log =
                              2))
)

## Print the results
svmFit # Based on the 10 cv's 3 times sigma = 0.009995619 C = 32 gives the most accuracy
svmFit$results
#plot average performance (i.e. ROC profile)
plot(svmFit, scales = list(x=list(log=2)))

## Test set predictions

svmPred <- predict(svmFit, finalTest, type="prob")
str(svmPred)

#Get confusion Matrix
svmPred<-predict(svmFit,finalTest[,names(finalTest)!="Class"])
confusionMatrix(svmPred, finalTest$Class)




## Use the "type" option to get class probabilities

predictedProbs <- predict(svmFit, newdata = finalTest, type = "prob")


#Print ROC curve
library(pROC)
library(ROCR)
svmPred <- predict(svmFit, finalTest, type="prob")
svmROC <- roc(finalTest$Class, svmPred[,"Severe"])
svmROC







## Fit the same model using different resampling methods. The main syntax change
## is the control object.
# 
# set.seed(1056)
# svmFit10CV <- train(Class ~ .,
#                     data = finalTrain,
#                     method = "svmRadial",
#                     tuneGrid = svmTuneGrid,
#                     trControl = trainControl(method = "cv", number = 10))
# svmFit10CV
# plot(svmFit10CV,scales = list(x=list(log=2)))
# 
# 
# set.seed(1056)
# svmFitLOO <- train(Class ~ .,
#                    data = finalTrain,
#                    method = "svmRadial",
#                    tuneGrid = svmTuneGrid,
#                    trControl = trainControl(method = "LOOCV"))
# svmFitLOO
# plot(svmFitLOO,scales = list(x=list(log=2)))
# 
# 
# set.seed(1056)
# svmFitLGO <- train(Class ~ .,
#                    data = finalTrain,
#                    method = "svmRadial",
#                    tuneGrid = svmTuneGrid,
#                    trControl = trainControl(method = "LGOCV", 
#                                             number = 50, 
#                                             p = .8))
# svmFitLGO 
# plot(svmFitLGO,scales = list(x=list(log=2)))
# 
# set.seed(1056)
# svmFitBoot <- train(Class ~ .,
#                     data = finalTrain,
#                     method = "svmRadial",
#                     tuneGrid = svmTuneGrid,
#                     trControl = trainControl(method = "boot", number = 50))
# svmFitBoot
# plot(svmFitBoot,scales = list(x=list(log=2)))
# set.seed(1056)
# svmFitBoot632 <- train(Class ~ .,
#                        data = finalTrain,
#                        method = "svmRadial",
#                        tuneGrid = svmTuneGrid,
#                        trControl = trainControl(method = "boot632", 
#                                                 number = 50))
# svmFitBoot632
# plot(svmFitBoot632,scales = list(x=list(log=2)))
# 
# #Use MARS less opaque
# 
# marsGrid <- data.frame(degree = 1, nprune = (1:10) * 3)
# 
# marsFit <- train(Class ~ .,
#                        data = finalTrain,
#                        method = "earth",
#                        tuneGrid = marsGrid,
#                        trControl = trainControl(method = "repeatedcv", 
#                                                 repeats = 5,
#                                                 classProbs = TRUE,
#                                                 selectionFunction = "tolerance"))
# marsFit
# plot(marsFit)

# 4.5 glm simple model#####
#install.packages("glmnet")
install.packages("arm")
library(arm)
library(glmnet)

set.seed(1056)
cvCtrl$sampling <- "smote"
glmGrid <- expand.grid(.alpha = (1:10) * 0.1, 
                       .lambda = (1:10) * 0.1)
glmProfile <- train(Class ~ .,
                    data = finalTrain,
                    tuneGrid = glmGrid,
                    metric = "Sens",
                    method = "glmnet",
                    trControl = cvCtrl)
# bglmProfile <- train(Class ~ .,
#                     data = finalTrain,
#                     method = "bayesglm",
#                     trControl = cvCtrl)
glmProfile
plot(varImp(glmProfile,scale=F))
# bglmProfile
glmPred <- predict(glmProfile, newdata = finalTest, type = "prob")
# bglmPred <- predict(bglmProfile, newdata = finalTest, type = "prob")

glmROC <- roc(finalTest$Class, glmPred[,"Severe"])
# bglmROC <- roc(finalTest$Class, bglmPred[,"Severe"])

#5. Results #####
#Get confusion Matrix's of all algorithms
glmPred<-predict(glmProfile,finalTest[,names(finalTest)!="Class"])
# bglmPred<-predict(bglmProfile,finalTest[,names(finalTest)!="Class"])
svmPred<-predict(svmFit,finalTest[,names(finalTest)!="Class"])

confusionMatrix(glmPred, finalTest$Class) #glm
# confusionMatrix(bglmPred, finalTest$Class) #bglm
confusionMatrix(svmPred, finalTest$Class) #svm
confusionMatrix(rfPred,finalTest$Class) #random forest
confusionMatrix(rpartPred2,finalTest$Class) #CART


#AUC's 
rpartROC
rfROC
svmROC
glmROC
# bglmROC

#Plot all

plot(rpartROC,type="S", xlim=c(1,0), legacy.axes = TRUE, asp = NA)
plot(rfROC,add = TRUE, type="S", col="blue")
plot(svmROC, add=TRUE, type="S", col="red")
plot(glmROC, add=TRUE, type="S", col="green")
legend("bottomright", # places a legend at the appropriate place \
       c('CART','RandomForest', 'SVM','GLM'), # puts text in the legend
       lty=c(1,1,1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5,2.5,2.5),col=c('black','blue','red','green')) # gives the legend lines the correct color and width




resamp <- resamples(list(SVM = svmFit, Logistic = glmProfile))
summary(resamp)

## These results are slightly different from those shown in the text.
## There are some differences in the train() function since the 
## original results were produced. This is due to a difference in
## predictions from the ksvm() function when class probs are requested
## and when they are not. See, for example, 
## https://stat.ethz.ch/pipermail/r-help/2013-November/363188.html

modelDifferences <- diff(resamp)
summary(modelDifferences) #statistically significant difference between logistic and svm SVM better

#Paired-t-test
modelDifferences$statistics$Accuracy



#Evaluate model performance
#Use random forest and quadratic discriminant models
install.packages("randomForest")
library(randomForest)
rfModel <- randomForest(Class ~ .,
                        data=finalTrain,
                        ntree = 2000 )

rfTestPred <- predict(rfModel, finalTest, type = "prob")
head(rfTestPred)
finalTest$RFprob <- rfTestPred[,"Nonsevere"]
finalTest$RFclass <- predict(rfModel, finalTest)
finalTest$lrclass <- predict(glmProfile, finalTest)
lrTestPred <- predict(glmProfile, finalTest, type="prob")
finalTest$lrprob <- lrTestPred[,"Nonsevere"]

#Sensitivity and specificity
#Nonsevere will be used as event and severe as nonevent
sensitivity(data=finalTest$RFclass,
            reference = finalTest$Class,
            positive = "Nonsevere")
specificity(data=finalTest$RFclass,
            reference = finalTest$Class,
            negative = "Severe")

#Sensitivity 0.375
#specificity 0.862

#Positivie and negative predicted values
posPredValue(data=finalTest$RFclass,
            reference = finalTest$Class,
            positive = "Nonsevere")
negPredValue(data=finalTest$RFclass,
            reference = finalTest$Class,
            negative = "Severe")
posPredValue(data=finalTest$RFclass,
            reference = finalTest$Class,
            positive = "Nonsevere",
            prevalence = 0.9)

# PPV 0.428 given that they are tested positive, 42.8% non-severe
# NPV 0.833 given that they are tested negative, 83.3% severe
#If change prevalence to 0.9, then PPV increases to 0.961

#Confusion Matrix
confusionMatrix(data = finalTest$RFclass,
                reference = finalTest$Class,
                positive = "Severe") #For RF model
confusionMatrix(data = finalTest$lrclass,
                reference = finalTest$Class,
                positive = "Severe") #For LR model
#Sensitivity 82.7, specificity 0.5

#ROC curve
library(pROC)
rocCurve <- roc(response=finalTest$Class,
                predictor = finalTest$RFprob,
                levels = rev(levels(finalTest$Class)))

auc(rocCurve) #0.7026
ci.auc(rocCurve) #0.4821,0.923
plot(rocCurve, legacy.axes = TRUE)

postResample(pred = finalTest$RFclass, obs = finalTest$Class)
postResample(pred = finalTest$lrclass, obs = finalTest$Class)
twoClassSummary(finalTest, lev = levels(finalTest$Class))

##### K Nearest Neighbor Modedl 2B


***2B1***: KNN model 

```{r , echo=FALSE, message=FALSE, warnings=FALSE}
nzv <- nearZeroVar(total, saveMetrics= TRUE)
nzv[nzv$nzv,][1:10,]

nzv <- nearZeroVar(total)
filteredtotal <- total[, -nzv]

set.seed(1234)
splitIndex <- createDataPartition(filteredtotal[,outcomeName], p = .6, list = FALSE, times = 1)
training  <- filteredtotal[ splitIndex,]
testing   <- filteredtotal[-splitIndex,]

trainX <- training[,names(training) != "COHORT2"]
preProcValues <- preProcess(x = trainX,method = c("center", "scale"))
preProcValues
prop.table(table(training$COHORT2)) * 100

set.seed(1234)
ctrl <- trainControl(method="cv",repeats = 3) #,classProbs=TRUE,summaryFunction = twoClassSummary)
knnFit <- train(COHORT2 ~ ., data = training, method = "knn", trControl = ctrl, preProcess = c("center","scale"), tuneLength = 20)
knnFit
plot(knnFit)
```

We see that the number of neighgbors near 13 gives the highest accuracy with regards to the model prediction ability.


***2B2***: KNN Contingency/Confusion Matrix 

```{r , echo=FALSE, message=FALSE, warnings=FALSE}
confusionMatrix(knnPredict, testing$COHORT2 )
mean(knnPredict == testing$COHORT2)
```

We see that the number of neighgbors near 13 gives the highest accuracy with regards to the model prediction ability.

***2B3***: KNN 2 class Summary  

```{r , echo=FALSE, message=FALSE, warnings=FALSE}
ctrl <- trainControl(method="cv",repeats = 3 ,classProbs=TRUE,summaryFunction = twoClassSummary)
knnFit <- train(COHORT2 ~ ., data = training, method = "knn", trControl = ctrl, preProcess = c("center","scale"), tuneLength = 20)
plot(knnFit, print.thres = 0.5, type="S")
```

The ROC for the cross-validated sets indicate that KNN isn't particularly good with recards to discriminating the different severity when the nearest neighbors bandwidth is too large.

***2B4***: KNN 2 ROC accuracy  

```{r , echo=FALSE, message=FALSE, warnings=FALSE}
confusionMatrix(knnPredict, testing$COHORT2 )
mean(knnPredict == testing$COHORT2)
```

The ROC for the cross-validated sets indicate that KNN isn't particularly good with recards to discriminating the different severity.


***2B4***: KNN 2 ROC accuracy  

```{r , echo=FALSE, message=FALSE, warnings=FALSE}
library(pROC)
knnPredict <- predict(knnFit,newdata = testing , type="prob")
knnROC <- roc(testing$Direction,knnPredict[,"Down"], levels = rev(testing$Direction))
knnROC
```

The ROC for the cross-validated sets indicate that KNN isn't particularly good with recards to discriminating the different severity.






***2c***: Evaluation of KNN model

```{r , echo=FALSE, message=FALSE, warnings=FALSE}
library(ROCR)

p <- predict(model, testDF, type="response")
pr <- prediction(p, testDF$COHORT3)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)


auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
```
The KNN model has an AUC of `r auc`.






##### Step 2: Model Training and Tuning

#####GMB Modedl 2A

***2a1***: Basic Parameter Tuning

```{r , echo=FALSE, message=FALSE, warnings=FALSE}
objControl  <- trainControl(method='cv', 
number=3, 
returnResamp='none', 
summaryFunction = twoClassSummary, 
classProbs = TRUE)
```

***2a2***: GBM training

```{r , echo=FALSE, message=FALSE, warnings=FALSE}
objModel <- train(trainDF[,predictorsNames], trainDF[,outcomeName], 
method='gbm', 
trControl=objControl,  
metric = "ROC",
preProc = c("center", "scale"))

summary(objModel)
print(objModel)
```


***2a3***: Evaluation of GBM model

```{r , echo=FALSE, message=FALSE, warnings=FALSE}
predictions <- predict(object=objModel, testDF[,predictorsNames], type='raw')
print(postResample(pred=predictions, obs=as.factor(testDF[,outcomeName])))
gbmac<-as.data.frame(postResample(pred=predictions, obs=as.factor(testDF[,outcomeName]))) 

```
Based on the results above, the GBM model is correct`r gbmac[1,1]`% of the time.   


```{r , echo=FALSE, message=FALSE, warnings=FALSE}
library(pROC)
library(ROCR)
predictions <- predict(object=objModel, testDF[,predictorsNames], type='prob')
head(predictions)

auc <- roc(ifelse(testDF[,outcomeName]=="Severe",1,0), predictions[[2]])
print(auc$auc)
gbmauc<-as.data.frame(auc$auc)
gbmauc[1,1]
plot(auc, main="ROC for GBM")
```

The GMB model has an AUC of `r gbmac[1,1]`.
