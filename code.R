###library
#package for simulating data
library(Umpire)
#packages for employed methods
library(RLT)
library(BART)
library(randomForest)
#package for application data
library(datamicroarray)
#package for cross validation
library(caret)
#package for parallel calculations
library(foreach)
library(doParallel)




##############################DATA SIMULATION FUNCION#####################################
###########################################################################################
########size of block - sizeb - {5, 15, 30}#################################################
#######number of observations - nobs - {50, 100, 200}#####################################



data_sim<-function(sizeb,nobs){
checkneg=0  #have to manually check whether in the train sample there are at least 20% positive cases and 20% negative cases
checkpos=0
while((checkpos<0.2*nobs)|(checkneg<0.2*nobs)){  
################################train data####################################
## Build a CancerModel with 1 subtype
nBlocks <- 20 # number of possible hits
cm <- CancerModel(name="cansim",
                  nPossible=nBlocks,
                  nPattern=1,
                  OUT = function(n) rnorm(n, 0, 1))
## Include 100 blocks/pathways that are not hit by cancer
nTotalBlocks <- nBlocks + 100
## Assign values to hyperparameters
## block size
blockSize <- round(rnorm(nTotalBlocks, sizeb, 0.3*sizeb))
blockSize<-replace(blockSize,blockSize<=0 ,1)
## log normal mean hypers
mu0 <- 6
sigma0 <- 1.5
## log normal sigma hypers
rate <- 28.11
shape <- 44.25
## block corr
p <- 0.6
w <- 5
#transcriptional activity
active <- 0.7
## Set up the baseline Engine
rho <- rbeta(nTotalBlocks, p*w, (1-p)*w)
base <- lapply(1:nTotalBlocks,
               function(i) {
                 bs <- blockSize[i]
                 co <- matrix(rho[i], nrow=bs, ncol=bs)
                 diag(co) <- 1
                 mu <- rnorm(bs, mu0, sigma0)
                 sigma <- matrix(1/rgamma(bs, rate=rate, shape=shape), nrow=1)
                 covo <- co *(t(sigma) %*% sigma)
                 MVN(mu, covo)
               })
eng <- EngineWithActivity(active, base, 2)
## Alter the means if there is a hit
altered <- alterMean(eng, normalOffset, delta=0, sigma=1)
## Build the CancerEngine using character strings
object <- CancerEngine(cm, "eng", "altered")
#summary(object)
## Simulate the data
temp <- rand(object, nobs+100)


## Add noise
nu <- 10
tau <- 20
phi <- 0.1
nm <- NoiseModel(nu, tau, phi)
temp$data<- blur(nm, temp$data) #realData


indx<-sample(1:(nobs+100),nobs)
dset<-list()
dset$data<-temp$data[,indx]
dset$clinical<-temp$clinical[indx,]


checkpos=sum(dset$clinical=="Bad")
checkneg=sum(dset$clinical=="Good")
}
###############################test data#######################################
#dsettest <- rand(object, 100)
#dsettest$data<- blur(nm, dsettest$data) #realData

dsettest<-list()
dsettest$data<-temp$data[,-indx]
dsettest$clinical<-temp$clinical[-indx,]



returnlist=list("dset"=dset,"dsettest"=dsettest)

return(returnlist)

}


##############################SIMULATION EXPERIMENT FUNCION################################
#100 repetitions for each setting of parameters for simulated data
###########################################################################################

simulation<-function(n,sb_k,reps){
#container for calculating TP,TN,FP, FN for each method
ERT=matrix(NA,ncol=5,nrow=reps,dimnames = list(c(1:reps),c('TP','TN','FP','FN','accur')))
BART=matrix(NA,ncol=5,nrow=reps,dimnames = list(c(1:reps),c('TP','TN','FP','FN','accur')))
DART=matrix(NA,ncol=5,nrow=reps,dimnames = list(c(1:reps),c('TP','TN','FP','FN','accur')))
RF=matrix(NA,ncol=5,nrow=reps,dimnames = list(c(1:reps),c('TP','TN','FP','FN','accur')))

  


for (i in 1:reps){

#simulating data with given parameters    
temp=data_sim(sizeb=sb_k,nobs=n)
temp_data = temp$dset
temp_datatest=temp$dsettest

#decoding
temp_data$clinical[,2]=ifelse(temp_data$clinical[,2]=='Bad',1,0)
temp_datatest$clinical[,2]=ifelse(temp_datatest$clinical[,2]=='Bad',1,0)

#train data
trainx = t(temp_data$data)
trainyfact = t(as.factor(as.matrix(t(temp_data$clinical[,2]))))
trainy = t(as.matrix(t(temp_data$clinical[,2])))

#test data
testx = t(temp_datatest$data)
testyfact= t(as.factor(as.matrix(t(temp_datatest$clinical[,2]))))
testy = t(as.matrix(t(temp_datatest$clinical[,2])))






#####################Random Forest######################################################


RF.fit <- randomForest(trainx,trainyfact,type='classification',mtry=dim(trainx)[2],nodesize=round(n^(1/3)))
RF.pred<-predict(RF.fit ,testx,type='response')


TP=sum((RF.pred=='1')&(testyfact=='1'))
TN=sum((RF.pred=='0')&(testyfact=='0'))
FP=sum((RF.pred=='1')&(testyfact=='0'))
FN=sum((RF.pred=='0')&(testyfact=='1'))


RF[i,'TP']=TP
RF[i,'TN']=TN
RF[i,'FP']=FP
RF[i,'FN']=FN
RF[i,'accur']=(TP+TN)/(TP+TN+FP+FN)



################Extremely Randomised Trees##########################################

ERT.fit = RLT(trainx, trainyfact, model = "classification")
ERT.pred = predict(ERT.fit, testx)

TP=sum((ERT.pred$Prediction=='1')&(testyfact=='1'))
TN=sum((ERT.pred$Prediction=='0')&(testyfact=='0'))
FP=sum((ERT.pred$Prediction=='1')&(testyfact=='0'))
FN=sum((ERT.pred$Prediction=='0')&(testyfact=='1'))


ERT[i,'TP']=TP
ERT[i,'TN']=TN
ERT[i,'FP']=FP
ERT[i,'FN']=FN
ERT[i,'accur']=(TP+TN)/(TP+TN+FP+FN)


#################################BART################################################


post<-lbart(trainx, trainy, nskip=100, ntree=20, ndpost=1000,sparse=FALSE,binaryOffset=0,k=3,power=10,base=0.75,rho=dim(trainx)[2],augment=TRUE)
pred <- predict(post, testx)
pred_y<-ifelse(pred$prob.test.mean<0.5,0,1)


TP=sum((pred_y==1)&(testy==1))
NT=sum((pred_y==0)&(testy==0))
FP=sum((pred_y==1)&(testy==0))
FN=sum((pred_y==0)&(testy==1))


BART[i,'TP']=TP
BART[i,'TN']=TN
BART[i,'FP']=FP
BART[i,'FN']=FN
BART[i,'accur']=(TP+TN)/(TP+TN+FP+FN)


#################################DART################################################


post<-lbart(trainx, trainy, nskip=100, ntree=50, ndpost=1000,sparse=TRUE,binaryOffset=0,k=3,power=10,base=0.75,a=1,rho=dim(trainx)[2],augment=TRUE)
pred <- predict(post, testx)
pred_y<-ifelse(pred$prob.test.mean<0.5,0,1)


TP=sum((pred_y==1)&(testy==1))
TN=sum((pred_y==0)&(testy==0))
FP=sum((pred_y==1)&(testy==0))
FN=sum((pred_y==0)&(testy==1))


DART[i,'TP']=TP
DART[i,'TN']=TN
DART[i,'FP']=FP
DART[i,'FN']=FN
DART[i,'accur']=(TP+TN)/(TP+TN+FP+FN)



}

#TP
TPBART=BART[,"TP"]
TPDART=DART[,"TP"]
TPERT=ERT[,"TP"]
TPRF=RF[,"TP"]
#TN
TNBART=BART[,"TN"]
TNDART=DART[,"TN"]
TNERT=ERT[,"TN"]
TNRF=RF[,"TN"]
#FP
FPBART=BART[,"FP"]
FPDART=DART[,"FP"]
FPERT=ERT[,"FP"]
FPRF=RF[,"FP"]

#FN
FNBART=BART[,"FN"]
FNDART=DART[,"FN"]
FNERT=ERT[,"FN"]
FNRF=RF[,"FN"]

return_list=list("TPBART"=TPBART,"TPDART"=TPDART,
                 "TPERT"=TPERT,"TPRF"=TPRF,
                 "TNBART"=TNBART,"TNDART"=TNDART,
                 "TNERT"=TNERT,"TNRF"=TNRF,
                 "FPBART"=FPBART,"FPDART"=FPDART,
                 "FPERT"=FPERT,"FPRF"=FPRF,
                 "FNBART"=FNBART,"FNDART"=FNDART,
                 "FNERT"=FNERT,"FNRF"=FNRF)


return(return_list)

}



############################IMPLEMENTATION OF SIMULATION EXPERIMENT########################
###########################################################################################
###########################################################################################


#5 - 600
#15 - 1500
#30 - 3500


list2=c("TPBART","TPDART","TPERT","TPRF",
        "TNBART","TNDART","TNERT","TNRF",
        "FPBART","FPDART","FPERT","FPRF",
        "FNBART","FNDART","FNERT","FNRF")


simulations_res5<-list()
simulations_res15<-list()
simulations_res30<-list()

#seting for parallel calculations
no_cores=detectCores()-1
cl <- makePSOCKcluster(no_cores, outfile='')
registerDoParallel(cl)
#registerDoSEQ()
#stopCluster(cl)
reps_n=100 #number of repetitions

Mark0 = proc.time()
simulations_res5 = foreach(n_i=c(50,100,200),.combine=rbind, .packages = c("dplyr",'Umpire','randomForest',"RLT","BART")) %dopar% {
set.seed(391629+n_i)
sb_k_i=5
simulations5=matrix(NA,ncol=16,nrow=reps_n,dimnames = list(c(1:reps_n),list2))
sim<-simulation(n=n_i,sb_k=sb_k_i,reps=reps_n)  

simulations5[,"TPBART"]<-sim$TPBART
simulations5[,"TPDART"]<-sim$TPDART
simulations5[,"TPERT"]<-sim$TPERT
simulations5[,"TPRF"]<-sim$TPRF

simulations5[,"TNBART"]<-sim$TNBART
simulations5[,"TNDART"]<-sim$TNDART
simulations5[,"TNERT"]<-sim$TNERT
simulations5[,"TNRF"]<-sim$TNRF

simulations5[,"FPBART"]<-sim$FPBART
simulations5[,"FPDART"]<-sim$FPDART
simulations5[,"FPERT"]<-sim$FPERT
simulations5[,"FPRF"]<-sim$FPRF

simulations5[,"FNBART"]<-sim$FNBART
simulations5[,"FNDART"]<-sim$FNDART
simulations5[,"FNERT"]<-sim$FNERT
simulations5[,"FNRF"]<-sim$FNRF

simulations5

}
#85 mins approx

simulations_res15 = foreach(n_i=c(50,100,200),.combine=rbind, .packages = c("dplyr",'Umpire','randomForest',"RLT","BART")) %dopar% {
  set.seed(391629+n_i)
  sb_k_i=15
  simulations15=matrix(NA,ncol=16,nrow=reps_n,dimnames = list(c(1:reps_n),list2))
  sim<-simulation(n=n_i,sb_k=sb_k_i,reps=100)  
  
  
  simulations15[,"TPBART"]<-sim$TPBART
  simulations15[,"TPDART"]<-sim$TPDART
  simulations15[,"TPERT"]<-sim$TPERT
  simulations15[,"TPRF"]<-sim$TPRF
  
  simulations15[,"TNBART"]<-sim$TNBART
  simulations15[,"TNDART"]<-sim$TNDART
  simulations15[,"TNERT"]<-sim$TNERT
  simulations15[,"TNRF"]<-sim$TNRF
  
  simulations15[,"FPBART"]<-sim$FPBART
  simulations15[,"FPDART"]<-sim$FPDART
  simulations15[,"FPERT"]<-sim$FPERT
  simulations15[,"FPRF"]<-sim$FPRF
  
  simulations15[,"FNBART"]<-sim$FNBART
  simulations15[,"FNDART"]<-sim$FNDART
  simulations15[,"FNERT"]<-sim$FNERT
  simulations15[,"FNRF"]<-sim$FNRF
  
  simulations15
  
}

#125 mins approx

simulations_res30 = foreach(n_i=c(50,100,200),.combine=rbind, .packages = c("dplyr",'Umpire','randomForest',"RLT","BART")) %dopar% {
  set.seed(391629+n_i)
  sb_k_i=30
  simulations30=matrix(NA,ncol=16,nrow=reps_n,dimnames = list(c(1:reps_n),list2))
  
  sim<-simulation(n=n_i,sb_k=sb_k_i,reps=100)  
  
  
  simulations30[,"TPBART"]<-sim$TPBART
  simulations30[,"TPDART"]<-sim$TPDART
  simulations30[,"TPERT"]<-sim$TPERT
  simulations30[,"TPRF"]<-sim$TPRF
  
  simulations30[,"TNBART"]<-sim$TNBART
  simulations30[,"TNDART"]<-sim$TNDART
  simulations30[,"TNERT"]<-sim$TNERT
  simulations30[,"TNRF"]<-sim$TNRF
  
  simulations30[,"FPBART"]<-sim$FPBART
  simulations30[,"FPDART"]<-sim$FPDART
  simulations30[,"FPERT"]<-sim$FPERT
  simulations30[,"FPRF"]<-sim$FPRF
  
  simulations30[,"FNBART"]<-sim$FNBART
  simulations30[,"FNDART"]<-sim$FNDART
  simulations30[,"FNERT"]<-sim$FNERT
  simulations30[,"FNRF"]<-sim$FNRF
  
  simulations30
  
}
#####210 mins approx

timeH20=proc.time() - Mark0
timeH20
#######################

#overall 432 mins (7 hours 12 mins) approx
#25932.97 secs
#3 cores were used (can be faster with more cores)



##########################EVALUATING THE PERFORMANCE#######################################
###########################################################################################
##########################################################################################


#function for caclulating the evaluation measures for each method
#precision, recall, F1, specificity
evaluation<-function(data){
  
#precision  
precRF=data[,'TPRF']/(data[,'TPRF']+data[,'FPRF'])
precERT=data[,'TPERT']/(data[,'TPERT']+data[,'FPERT'])
precBART=data[,'TPBART']/(data[,'TPBART']+data[,'FPBART'])
precDART=data[,'TPDART']/(data[,'TPDART']+data[,'FPDART'])
#recall
recRF=data[,'TPRF']/(data[,'TPRF']+data[,'FNRF'])
recERT=data[,'TPERT']/(data[,'TPERT']+data[,'FNERT'])
recBART=data[,'TPBART']/(data[,'TPBART']+data[,'FNBART'])
recDART=data[,'TPDART']/(data[,'TPDART']+data[,'FNDART'])
#F1
F1RF=2*precRF*recRF/(precRF+recRF)
F1ERT=2*precERT*recERT/(precERT+recERT)
F1BART=2*precBART*recBART/(precBART+recBART)
F1DART=2*precDART*recDART/(precDART+recDART)

#accuracy
accRF=(data[,'TPRF']+data[,'TNRF'])/(data[,'FPRF']+data[,'FNRF']+data[,'TPRF']+data[,'TNRF'])
accERT=(data[,'TPERT']+data[,'TNERT'])/(data[,'FPERT']+data[,'FNERT']+data[,'TPERT']+data[,'TNERT'])
accBART=(data[,'TPBART']+data[,'TNBART'])/(data[,'FPBART']+data[,'FNBART']+data[,'TPBART']+data[,'TNBART'])
accDART=(data[,'TPDART']+data[,'TNDART'])/(data[,'FPDART']+data[,'FNDART']+data[,'TPDART']+data[,'TNDART'])

return_list=list('recRF'=recRF,'recERT'=recERT,'recBART'=recBART,'recDART'=recDART,
'precRF'=precRF,'precERT'=precERT,'precBART'=precBART,'precDART'=precDART,
'F1RF'=F1RF,'F1ERT'=F1ERT,'F1BART'=F1BART,'F1DART'=F1DART,
'accRF'=accRF,'accERT'=accERT,'accBART'=accBART,'accDART'=accDART)

return(return_list)
}




ev5<-evaluation(simulations_res5)
ev15<-evaluation(simulations_res15)
ev30<-evaluation(simulations_res30)

###################APPLICATION FOR A REAL DATA#################################
###############################################################################
###############################################################################

#downloading real data
data('gravier', package = 'datamicroarray')
gravier$y<-as.matrix(gravier$y)



#evaluate performance with 5-fold cross validation
folds=10

ERT=matrix(NA,ncol=5,nrow=folds,dimnames = list(c(1:folds),c('TP','TN','FP','FN','accur')))
BART=matrix(NA,ncol=5,nrow=folds,dimnames = list(c(1:folds),c('TP','TN','FP','FN','accur')))
DART=matrix(NA,ncol=5,nrow=folds,dimnames = list(c(1:folds),c('TP','TN','FP','FN','accur')))
RF=matrix(NA,ncol=5,nrow=folds,dimnames = list(c(1:folds),c('TP','TN','FP','FN','accur')))
  
#create folds
set.seed(391629)
flds <- createFolds(gravier$y, k = folds, list = TRUE, returnTrain = FALSE)


###running CV
CV_res=list()
CV_res=foreach(i=c(1:folds),.combine=rbind, .packages = c("dplyr",'Umpire','randomForest',"RLT","BART")) %dopar% {
names(flds)[i] <- "test"  
set.seed(391629+i)
  
#test data
testy<-as.matrix(gravier$y[flds$test,])
testy<-ifelse(testy=='poor',1,0)
testyfact<-as.factor(testy)
testx<-gravier$x[flds$test, ]

#train data
trainy<-as.matrix(gravier$y[-flds$test,])
trainy<-ifelse(trainy=='poor',1,0)
trainyfact<-as.factor(trainy)
trainx<-gravier$x[-flds$test, ]

###Random Forest
RF.fit <- randomForest(trainx,trainyfact,type='classification',mtry=dim(trainx)[2],nodesize=round((dim(trainx)[1])^(1/3)))
RF.pred<-predict(RF.fit ,testx,type='response')



TP=sum((RF.pred=='1')&(testyfact=='1'))
TN=sum((RF.pred=='0')&(testyfact=='0'))
FP=sum((RF.pred=='1')&(testyfact=='0'))
FN=sum((RF.pred=='0')&(testyfact=='1'))


RF[i,'TP']=TP
RF[i,'TN']=TN
RF[i,'FP']=FP
RF[i,'FN']=FN
RF[i,'accur']=(TP+TN)/(TP+TN+FP+FN)




###Extremely Randomized Trees
ERT.fit = RLT(trainx, trainyfact, model = "classification")
ERT.pred = predict(ERT.fit, testx)

TP=sum((ERT.pred$Prediction=='1')&(testyfact=='1'))
TN=sum((ERT.pred$Prediction=='0')&(testyfact=='0'))
FP=sum((ERT.pred$Prediction=='1')&(testyfact=='0'))
FN=sum((ERT.pred$Prediction=='0')&(testyfact=='1'))


ERT[i,'TP']=TP
ERT[i,'TN']=TN
ERT[i,'FP']=FP
ERT[i,'FN']=FN
ERT[i,'accur']=(TP+TN)/(TP+TN+FP+FN)


###BART
post<-lbart(trainx, trainy, nskip=100, ntree=20, ndpost=1000,sparse=FALSE,binaryOffset=0,k=3,power=3,base=0.9,rho=dim(trainx)[2],augment=TRUE)
pred <- predict(post, testx)
pred_y<-ifelse(pred$prob.test.mean<0.5,0,1)


TP=sum((pred_y==1)&(testy==1))
NT=sum((pred_y==0)&(testy==0))
FP=sum((pred_y==1)&(testy==0))
FN=sum((pred_y==0)&(testy==1))


BART[i,'TP']=TP
BART[i,'TN']=TN
BART[i,'FP']=FP
BART[i,'FN']=FN
BART[i,'accur']=(TP+TN)/(TP+TN+FP+FN)



###DART
post<-lbart(trainx, trainy, nskip=100, ntree=50, ndpost=1000,sparse=TRUE,binaryOffset=0,k=3,power=3,base=0.9,a=1,rho=dim(trainx)[2],augment=TRUE)
pred <- predict(post, testx)
pred_y<-ifelse(pred$prob.test.mean<0.5,0,1)


TP=sum((pred_y==1)&(testy==1))
TN=sum((pred_y==0)&(testy==0))
FP=sum((pred_y==1)&(testy==0))
FN=sum((pred_y==0)&(testy==1))


DART[i,'TP']=TP
DART[i,'TN']=TN
DART[i,'FP']=FP
DART[i,'FN']=FN
DART[i,'accur']=(TP+TN)/(TP+TN+FP+FN)


CV_temp<-cbind(BART[i,'TP'],DART[i,'TP'],ERT[i,'TP'],RF[i,'TP'],
              BART[i,'TN'],DART[i,'TN'],ERT[i,'TN'],RF[i,'TN'],
              BART[i,'FP'],DART[i,'FP'],ERT[i,'FP'],RF[i,'FP'],
              BART[i,'FN'],DART[i,'FN'],ERT[i,'FN'],RF[i,'FN'])

CV_temp
}
##approx 7 mins

#stop parallel execution
registerDoSEQ()
stopCluster(cl)



### evaluating
colnames(CV_res)=list2




#function
evaluation_CV<-function(data){
  #precision  
  precRF=mean(data$precRF,na.rm = TRUE)
  precERT=mean(data$precERT,na.rm = TRUE)
  precBART=mean(data$precBART,na.rm = TRUE)
  precDART=mean(data$precDART,na.rm = TRUE)
  #recall
  recRF=mean(data$recRF,na.rm = TRUE)
  recERT=mean(data$recERT,na.rm = TRUE)
  recBART=mean(data$recBART,na.rm = TRUE)
  recDART=mean(data$recDART,na.rm = TRUE)
  #F1
  F1RF=mean(data$F1RF,na.rm = TRUE)
  F1ERT=mean(data$F1ERT,na.rm = TRUE)
  F1BART=mean(data$F1BART,na.rm = TRUE)
  F1DART=mean(data$F1DART,na.rm = TRUE)
  #accuracy
  accRF=mean(data$accRF,na.rm = TRUE)
  accERT=mean(data$accERT,na.rm = TRUE)
  accBART=mean(data$accBART,na.rm = TRUE)
  accDART=mean(data$accDART,na.rm = TRUE)
  

  return_list=list('recRF'=round(recRF,4),'recERT'=round(recERT,4),'recBART'=round(recBART,4),'recDART'=round(recDART,4),
                   'precRF'=round(precRF,4),'precERT'=round(precERT,4),'precBART'=round(precBART,4),'precDART'=round(precDART,4),
                   'F1RF'=round(F1RF,4),'F1ERT'=round(F1ERT,4),'F1BART'=round(F1BART,4),'F1DART'=round(F1DART,4),
                   'accRF'=round(accRF,4),'accERT'=round(accERT,4),'accBART'=round(accBART,4),'accDART'=round(accDART,4))
  
  return(return_list)
}



###########################FORMATTING THE OUTPUT (TABLES)#############
#################################################################################
#################################################################################



####for Application

CV_ev<-evaluation_CV(evaluation(CV_res))


###organising output
colrec<-rbind(CV_ev$recBART,CV_ev$recDART,CV_ev$recRF)
colprec<-rbind(CV_ev$precBART,CV_ev$precDART,CV_ev$precRF)
colf1<-rbind(CV_ev$F1BART,CV_ev$F1DART,CV_ev$F1RF)
colacc<-rbind(CV_ev$accBART,CV_ev$accDART,CV_ev$accRF)
#colspec<-rbind(CV_ev$specBART,CV_ev$specDART,CV_ev$specRF)
CV_results<-cbind(colrec,colprec,colf1,colacc)
colnames(CV_results)<-c('recall','precision','F1','accuracy')
rownames(CV_results)<-c('BART','DART','RF')

for (i in 1:dim(CV_results)[2]){
  CV_results[,i] = replace(CV_results[,i],CV_results[,i]==max(CV_results[,i]),paste("!",as.character(max(CV_results[,i])),"!"))
}

#how to provide output
#CV_results



####For Simulations
#function for producing tables
result_table<-function(data){
  temp_table50=evaluation_CV(evaluation(data[1:100,]))
  temp_table100=evaluation_CV(evaluation(data[101:200,]))
  temp_table200=evaluation_CV(evaluation(data[201:300,]))
  
  col50<-cbind(rbind(temp_table50$recBART,temp_table50$recDART,temp_table50$recRF),
                  rbind(temp_table50$precBART,temp_table50$precDART,temp_table50$precRF),
                  rbind(temp_table50$F1BART,temp_table50$F1DART,temp_table50$F1RF),
               rbind(temp_table50$accBART,temp_table50$accDART,temp_table50$accRF))
  col100<-cbind(rbind(temp_table100$recBART,temp_table100$recDART,temp_table100$recRF),
               rbind(temp_table100$precBART,temp_table100$precDART,temp_table100$precRF),
               rbind(temp_table100$F1BART,temp_table100$F1DART,temp_table100$F1RF),
               rbind(temp_table100$accBART,temp_table100$accDART,temp_table100$accRF))
  col200<-cbind(rbind(temp_table200$recBART,temp_table200$recDART,temp_table200$recRF),
               rbind(temp_table200$precBART,temp_table200$precDART,temp_table200$precRF),
               rbind(temp_table200$F1BART,temp_table200$F1DART,temp_table200$F1RF),
               rbind(temp_table200$accBART,temp_table200$accDART,temp_table200$accRF))
  
  
  temp_results<-rbind(c('recall','precision','F1','accuracy',
                        'recall','precision','F1','accuracy',
                        'recall','precision','F1','accuracy'),
                      cbind(col50,col100,col200))
  #adding "!" around the best performance
  for (i in 1:dim(temp_results)[2]){
    temp_results[,i] = replace(temp_results[,i],temp_results[,i]==max(temp_results[2:4,i]),paste("!",as.character(max(temp_results[2:4,i])),"!"))
  }
  
  colnames(temp_results)<-c('N=50','','','','N=100','','','','N=200','','','')
  rownames(temp_results)<-c('','BART','DART','RF')
  
  
  return(temp_results)    
}


#how to use
#result_table(simulations_res15)

