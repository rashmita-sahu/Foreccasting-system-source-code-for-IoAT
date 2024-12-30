# Load necessary libraries
library(openssl)
library(DChaos)
library(Metrics)
library(cluster)
library(factoextra)

# Load additional libraries from your original code
library(VineCopula)
library(moments)
library(copula)
library(fCopulae)
library(ggplot2)
library(writexl)
library(sets)
library(mclust)
library(entropy)
library(Matrix)
library(FuzzyR)
library(GLDEX)
install.packages("genlasso")
library(genlasso)
# ==============================================
# PART 1: Load and preprocess the dataset
options(scipen=999)
Sys.time()
x<-read.csv('D:/R_Experiments/Irrigation_DS.csv')

x=c(x[ ,1], x[ ,2], x[ ,3])
y=x[ ,1]
z=y[1:50000]
############# Total variation Regularization on dataset##################################
xsignal = c(z)
output = admm.tv(c(z))
output=c(output)
## visualize
opar <- par(no.readonly=TRUE)
plot(xsignal, type="l", main="Total Variance Regularization for soil moisture Sensor")
lines(1:50000, output$x, col="red", lwd=1)
par(opar)
#####################################Reading Dataset for probablistic clustering##########################################################
x1=x[ , 1] #soil moisture sensor data
x2=x[ , 2] #temperature sensor data
x3=x[ , 3] #humidity sensor data

P_cluster<-Mclust(unlist(prop.table(x1)))
summary(P_cluster)
results = data.frame(unlist(prop.table(x1)),cluster=P_cluster$classification)
clusplot(as.matrix(data.frame(unlist(prop.table(x1)))), results$cluster, shade = FALSE,labels=5,col.clus="blue",col.p="red",span=FALSE,main="Probabilistic Cluster Mapping",cex=1.2)

v <- results$cluster
c1=0
c2=0
c3=0
c4=0
c5=0
c6=0
c7=0
c8=0
c9=0
for(i in 1:length(x1))
{
  
  if (v[i] == 1)
  {
    c1[i] <- x1[i]
    
  }
  if (v[i] == 2)
  {
    c2[i] <- x1[i]
    
  }
  if (v[i] == 3)
  {
    c3[i] <- x1[i]
    
  }
  if (v[i] == 4)
  {
    c4[i] <- x1[i]
    
  }
  if (v[i] == 5)
  {
    c5[i] <- x1[i]
    
  }
  if (v[i] == 6)
  {
    c6[i] <- x1[i]
    
  }
  if (v[i] == 7)
  {
    c7[i] <- x1[i]
    
  }
  if (v[i] == 8)
  {
    c8[i] <- x1[i]
    
  }
  if (v[i] == 9)
  {
    c9[i] <- x1[i]
    
  }
}

# Create a loop to iterate over indices 1 to 9
for (i in 1:9) {
  var_name <- paste0("c", i)  # Dynamically create variable name "c1", "c2", ..., "c9"
  
  # Use get() to access the variable by name, apply na.omit() and fun.zero.omit(), then assign it back
  assign(var_name, fun.zero.omit(as.vector(na.omit(get(var_name)))))
  assign(var_name, fun.zero.omit(as.vector(na.omit(get(var_name)))))
  
}

#######################Finding significant value from each cluster##############
i_value_c1=-log(prop.table(c1))
i_value_c2=-log(prop.table(c2))
i_value_c3=-log(prop.table(c3))
i_value_c4=-log(prop.table(c4))
i_value_c5=-log(prop.table(c5))
i_value_c6=-log(prop.table(c6))
i_value_c7=-log(prop.table(c7))
i_value_c8=-log(prop.table(c8))
i_value_c9=-log(prop.table(c9))

for(i in 1:length(c1))
{
  if(min(i_value_c1)==i_value_c1[i])
    sensitive_value_c1=c1[i]
}

for(i in 1:length(c2))
{
  if(min(i_value_c2)==i_value_c2[i])
    sensitive_value_c2=c2[i]
}

for(i in 1:length(c3))
{
  if(min(i_value_c3)==i_value_c3[i])
    sensitive_value_c3=c3[i]
}

for(i in 1:length(c4))
{
  if(min(i_value_c4)==i_value_c4[i])
    sensitive_value_c4=c4[i]
}

for(i in 1:length(c5))
{
  if(min(i_value_c5)==i_value_c5[i])
    sensitive_value_c5=c5[i]
}

for(i in 1:length(c6))
{
  if(min(i_value_c6)==i_value_c6[i])
    sensitive_value_c6=c6[i]
}

for(i in 1:length(c7))
{
  if(min(i_value_c7)==i_value_c7[i])
    sensitive_value_c7=c7[i]
}

for(i in 1:length(c8))
{
  if(min(i_value_c8)==i_value_c8[i])
    sensitive_value_c8=c8[i]
}

for(i in 1:length(c9))
{
  if(min(i_value_c9)==i_value_c9[i])
    sensitive_value_c9=c9[i]
}

# Loop through c1 to c9
for (i in 1:9) {
  var_name <- paste0("c", i)  # Dynamically create variable names "c1", "c2", ..., "c9"
  
  # Calculate i_value for c1 to c9
  i_value <- -log(prop.table(get(var_name)))
  assign(paste0("i_value_", var_name), i_value)
  
  # Find the sensitive value for each c1 to c9
  sensitive_value <- get(var_name)[which.min(i_value)]
  assign(paste0("sensitive_value_", var_name), sensitive_value)
}

#########################Lagrangian point Redefining##################################

R=(sensitive_value_c1-sensitive_value_c2)^2
r=R*(((sensitive_value_c1)/(3*sensitive_value_c2))^1/3)
conscious_value_for_c1_and_c2=sensitive_value_c1+r
conscious_value_for_c1_and_c2

R=(sensitive_value_c1-sensitive_value_c3)^2
r=R*(((sensitive_value_c1)/(3*sensitive_value_c3))^1/3)
conscious_value_for_c1_and_c3=sensitive_value_c1+r
conscious_value_for_c1_and_c3

R=(sensitive_value_c1-sensitive_value_c4)^2
r=R*(((sensitive_value_c1)/(3*sensitive_value_c4))^1/3)
conscious_value_for_c1_and_c4=sensitive_value_c1+r
conscious_value_for_c1_and_c4

R=(sensitive_value_c8-sensitive_value_c9)^2
r=R*(((sensitive_value_c8)/(3*sensitive_value_c9))^1/3)
conscious_value_for_c8_and_c9=sensitive_value_c8+r
conscious_value_for_c8_and_c9

############################Soil Temperature sensor values#############################
P_cluster1<-Mclust(unlist(prop.table(x2)))

#P_cluster<-Mclust(unlist(x1))
#summary(P_cluster)
#plot(P_cluster)
summary(P_cluster1)
#plot(P_cluster)
results1 = data.frame(unlist(prop.table(x2)),cluster1=P_cluster1$classification)
clusplot(as.matrix(data.frame(unlist(prop.table(x2)))), results1$cluster1, shade = FALSE,labels=5,col.clus="blue",col.p="red",span=FALSE,main="Probabilistic Cluster Mapping",cex=1.2)

v <- results1$cluster1
c1=0
c2=0
c3=0
c4=0
c5=0
c6=0
c7=0
c8=0
c9=0
for(i in 1:length(x2))
{
  
  if (v[i] == 1)
  {
    c1[i] <- x2[i]
    
  }
  if (v[i] == 2)
  {
    c2[i] <- x2[i]
    
  }
  if (v[i] == 3)
  {
    c3[i] <- x2[i]
    
  }
  if (v[i] == 4)
  {
    c4[i] <- x2[i]
    
  }
  if (v[i] == 5)
  {
    c5[i] <- x2[i]
    
  }
  if (v[i] == 6)
  {
    c6[i] <- x2[i]
    
  }
  if (v[i] == 7)
  {
    c7[i] <- x2[i]
    
  }
  if (v[i] == 8)
  {
    c8[i] <- x2[i]
    
  }
  if (v[i] == 9)
  {
    c9[i] <- x2[i]
    
  }
}

#######################Finding significant value from each cluster##############
i_value_c1=-log(prop.table(c1))
i_value_c2=-log(prop.table(c2))
i_value_c3=-log(prop.table(c3))
i_value_c4=-log(prop.table(c4))
i_value_c5=-log(prop.table(c5))
i_value_c6=-log(prop.table(c6))
i_value_c7=-log(prop.table(c7))
i_value_c8=-log(prop.table(c8))
i_value_c9=-log(prop.table(c9))

for(i in 1:length(c1))
{
  if(min(i_value_c1)==i_value_c1[i])
    sensitive_value_c1=c1[i]
}

for(i in 1:length(c2))
{
  if(min(i_value_c2)==i_value_c2[i])
    sensitive_value_c2=c2[i]
}

for(i in 1:length(c3))
{
  if(min(i_value_c3)==i_value_c3[i])
    sensitive_value_c3=c3[i]
}

for(i in 1:length(c4))
{
  if(min(i_value_c4)==i_value_c4[i])
    sensitive_value_c4=c4[i]
}

for(i in 1:length(c5))
{
  if(min(i_value_c5)==i_value_c5[i])
    sensitive_value_c5=c5[i]
}

for(i in 1:length(c6))
{
  if(min(i_value_c6)==i_value_c6[i])
    sensitive_value_c6=c6[i]
}

for(i in 1:length(c7))
{
  if(min(i_value_c7)==i_value_c7[i])
    sensitive_value_c7=c7[i]
}

for(i in 1:length(c8))
{
  if(min(i_value_c8)==i_value_c8[i])
    sensitive_value_c8=c8[i]
}

for(i in 1:length(c9))
{
  if(min(i_value_c9)==i_value_c9[i])
    sensitive_value_c9=c9[i]
}

#########################Lagrangian point Redefining##################################

R=(sensitive_value_c8-sensitive_value_c9)^2
r=R*(((sensitive_value_c8)/(3*sensitive_value_c9))^1/3)
conscious_value_for_c8_and_c9=sensitive_value_c8+r
conscious_value_for_c8_and_c9
############################Soil Humidity sensor Values#########################################

P_cluster2<-Mclust(unlist(prop.table(x3)))

#P_cluster<-Mclust(unlist(x1))
#summary(P_cluster)
#plot(P_cluster)
summary(P_cluster2)
#plot(P_cluster)
results2 = data.frame(unlist(prop.table(x3)),cluster2=P_cluster2$classification)
clusplot(as.matrix(data.frame(unlist(prop.table(x3)))), results2$cluster2, shade = FALSE,labels=5,col.clus="blue",col.p="red",span=FALSE,main="Probabilistic Cluster Mapping",cex=1.2)

v <- results2$cluster2
c1=0
c2=0
c3=0
c4=0
c5=0
c6=0
c7=0
c8=0
c9=0
for(i in 1:length(x3))
{
  
  if (v[i] == 1)
  {
    c1[i] <- x3[i]
    
  }
  if (v[i] == 2)
  {
    c2[i] <- x3[i]
    
  }
  if (v[i] == 3)
  {
    c3[i] <- x3[i]
    
  }
  if (v[i] == 4)
  {
    c4[i] <- x3[i]
    
  }
  if (v[i] == 5)
  {
    c5[i] <- x3[i]
    
  }
  if (v[i] == 6)
  {
    c6[i] <- x3[i]
    
  }
  if (v[i] == 7)
  {
    c7[i] <- x3[i]
    
  }
  if (v[i] == 8)
  {
    c8[i] <- x3[i]
    
  }
  if (v[i] == 9)
  {
    c9[i] <- x3[i]
    
  }
}


#######################Finding significant value from each cluster##############
i_value_c1=-log(prop.table(c1))
i_value_c2=-log(prop.table(c2))
i_value_c3=-log(prop.table(c3))
i_value_c4=-log(prop.table(c4))
i_value_c5=-log(prop.table(c5))
i_value_c6=-log(prop.table(c6))
i_value_c7=-log(prop.table(c7))
i_value_c8=-log(prop.table(c8))
i_value_c9=-log(prop.table(c9))

for(i in 1:length(c1))
{
  if(min(i_value_c1)==i_value_c1[i])
    sensitive_value_c1=c1[i]
}

for(i in 1:length(c2))
{
  if(min(i_value_c2)==i_value_c2[i])
    sensitive_value_c2=c2[i]
}

for(i in 1:length(c3))
{
  if(min(i_value_c3)==i_value_c3[i])
    sensitive_value_c3=c3[i]
}

for(i in 1:length(c4))
{
  if(min(i_value_c4)==i_value_c4[i])
    sensitive_value_c4=c4[i]
}

for(i in 1:length(c5))
{
  if(min(i_value_c5)==i_value_c5[i])
    sensitive_value_c5=c5[i]
}

for(i in 1:length(c6))
{
  if(min(i_value_c6)==i_value_c6[i])
    sensitive_value_c6=c6[i]
}

for(i in 1:length(c7))
{
  if(min(i_value_c7)==i_value_c7[i])
    sensitive_value_c7=c7[i]
}

for(i in 1:length(c8))
{
  if(min(i_value_c8)==i_value_c8[i])
    sensitive_value_c8=c8[i]
}

for(i in 1:length(c9))
{
  if(min(i_value_c9)==i_value_c9[i])
    sensitive_value_c9=c9[i]
}

#########################Lagrangian point Redefining##################################

R=(sensitive_value_c8-sensitive_value_c9)^2
r=R*(((sensitive_value_c8)/(3*sensitive_value_c9))^1/3)
conscious_value_for_c8_and_c9=sensitive_value_c8+r
conscious_value_for_c8_and_c9
###################Plausoble value extraction##########################################
X <- read.csv('E:/Paper_5/Exp_Result_5.csv')
P=prop.table(prop.table(X[, 1]))
Q=prop.table(prop.table(X[, 2]))
R=prop.table(prop.table(X[, 3]))
Concious_pattern=expand.grid(P, Q, R)
write_xlsx(data.frame(Concious_pattern), "E:/Paper_5/plausible value.xlsx")
e=0
for(i in 1: nrow(Concious_pattern))
{
  e[i] = (unlist(Concious_pattern[i,1]))*(unlist(Concious_pattern[i,2]))*(unlist(Concious_pattern[i,2]))
  
}
write_xlsx(data.frame(e), "E:/Paper_5/plusible value.xlsx")

max(e)
for(i in 1:length(e))
{
  if(e[i]==max(e))
    j=i
}
Concious_pattern[j, ]
W=X[, 2]
A=prop.table((W))
for (i in 1: length(X[, 1]))
{
  if(0.05132291==A[i])
    k=i
}

Metrics::mape(mean(c1), 26.17241)
Metrics::mape(mean(c2), 26.17241)

########################Soil Moisture cross validation##################################

Metrics::mape(mean(c1), conscious_value_for_c1_and_c2)
Metrics::mape(mean(c2), conscious_value_for_c1_and_c2)
Metrics::smape(mean(c1), conscious_value_for_c1_and_c2)
Metrics::smape(mean(c2), conscious_value_for_c1_and_c2)
fuzzyr.accuracy(mean(c1), conscious_value_for_c1_and_c2)
fuzzyr.accuracy(mean(c2), conscious_value_for_c1_and_c2)

Metrics::mape(mean(c2), conscious_value_for_c2_and_c3)
Metrics::mape(mean(c3), conscious_value_for_c2_and_c3)
Metrics::smape(mean(c2), conscious_value_for_c2_and_c3)
Metrics::smape(mean(c3), conscious_value_for_c2_and_c3)
fuzzyr.accuracy(mean(c2), conscious_value_for_c2_and_c3)
fuzzyr.accuracy(mean(c3), conscious_value_for_c2_and_c3)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c4)
Metrics::mape(mean(c4), conscious_value_for_c3_and_c4)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c4)
Metrics::smape(mean(c4), conscious_value_for_c3_and_c4)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c4)
fuzzyr.accuracy(mean(c4), conscious_value_for_c3_and_c4)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c5)
Metrics::mape(mean(c5), conscious_value_for_c3_and_c5)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c5)
Metrics::smape(mean(c5), conscious_value_for_c3_and_c5)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c5)
fuzzyr.accuracy(mean(c5), conscious_value_for_c3_and_c5)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c6)
Metrics::mape(mean(c6), conscious_value_for_c3_and_c6)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c6)
Metrics::smape(mean(c6), conscious_value_for_c3_and_c6)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c6)
fuzzyr.accuracy(mean(c6), conscious_value_for_c3_and_c6)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c3_and_c7)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c3_and_c7)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c3_and_c7)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c3_and_c8)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c3_and_c8)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c3_and_c8)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c3_and_c9)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c3_and_c9)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c3_and_c9)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c5)
Metrics::mape(mean(c5), conscious_value_for_c4_and_c5)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c5)
Metrics::smape(mean(c5), conscious_value_for_c4_and_c5)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c5)
fuzzyr.accuracy(mean(c5), conscious_value_for_c4_and_c5)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c6)
Metrics::mape(mean(c6), conscious_value_for_c4_and_c6)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c6)
Metrics::smape(mean(c6), conscious_value_for_c4_and_c6)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c6)
fuzzyr.accuracy(mean(c6), conscious_value_for_c4_and_c6)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c4_and_c7)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c4_and_c7)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c4_and_c7)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c4_and_c8)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c4_and_c8)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c4_and_c8)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c4_and_c9)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c4_and_c9)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c4_and_c9)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c6)
Metrics::mape(mean(c6), conscious_value_for_c5_and_c6)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c6)
Metrics::smape(mean(c6), conscious_value_for_c5_and_c6)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c6)
fuzzyr.accuracy(mean(c6), conscious_value_for_c5_and_c6)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c5_and_c7)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c5_and_c7)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c5_and_c7)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c5_and_c8)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c5_and_c8)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c5_and_c8)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c5_and_c9)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c5_and_c9)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c5_and_c9)

Metrics::mape(mean(c6), conscious_value_for_c6_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c6_and_c7)
Metrics::smape(mean(c6), conscious_value_for_c6_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c6_and_c7)
fuzzyr.accuracy(mean(c6), conscious_value_for_c6_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c6_and_c7)

Metrics::mape(mean(c6), conscious_value_for_c6_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c6_and_c8)
Metrics::smape(mean(c6), conscious_value_for_c6_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c6_and_c8)
fuzzyr.accuracy(mean(c6), conscious_value_for_c6_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c6_and_c8)

Metrics::mape(mean(c6), conscious_value_for_c6_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c6_and_c9)
Metrics::smape(mean(c6), conscious_value_for_c6_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c6_and_c9)
fuzzyr.accuracy(mean(c6), conscious_value_for_c6_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c6_and_c9)

Metrics::mape(mean(c7), conscious_value_for_c7_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c7_and_c8)
Metrics::smape(mean(c7), conscious_value_for_c7_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c7_and_c8)
fuzzyr.accuracy(mean(c7), conscious_value_for_c7_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c7_and_c8)

Metrics::mape(mean(c7), conscious_value_for_c7_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c7_and_c9)
Metrics::smape(mean(c7), conscious_value_for_c7_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c7_and_c9)
fuzzyr.accuracy(mean(c7), conscious_value_for_c7_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c7_and_c9)

Metrics::mape(mean(c8), conscious_value_for_c8_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c8_and_c9)
Metrics::smape(mean(c8), conscious_value_for_c8_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c8_and_c9)
fuzzyr.accuracy(mean(c8), conscious_value_for_c8_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c8_and_c9)

############################Temperature Validation Matrices########################################

Metrics::mape(mean(c1), conscious_value_for_c1_and_c2)
Metrics::mape(mean(c2), conscious_value_for_c1_and_c2)
Metrics::smape(mean(c1), conscious_value_for_c1_and_c2)
Metrics::smape(mean(c2), conscious_value_for_c1_and_c2)
fuzzyr.accuracy(mean(c1), conscious_value_for_c1_and_c2)
fuzzyr.accuracy(mean(c2), conscious_value_for_c1_and_c2)

Metrics::mape(mean(c2), conscious_value_for_c2_and_c3)
Metrics::mape(mean(c3), conscious_value_for_c2_and_c3)
Metrics::smape(mean(c2), conscious_value_for_c2_and_c3)
Metrics::smape(mean(c3), conscious_value_for_c2_and_c3)
fuzzyr.accuracy(mean(c2), conscious_value_for_c2_and_c3)
fuzzyr.accuracy(mean(c3), conscious_value_for_c2_and_c3)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c4)
Metrics::mape(mean(c4), conscious_value_for_c3_and_c4)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c4)
Metrics::smape(mean(c4), conscious_value_for_c3_and_c4)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c4)
fuzzyr.accuracy(mean(c4), conscious_value_for_c3_and_c4)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c5)
Metrics::mape(mean(c5), conscious_value_for_c3_and_c5)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c5)
Metrics::smape(mean(c5), conscious_value_for_c3_and_c5)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c5)
fuzzyr.accuracy(mean(c5), conscious_value_for_c3_and_c5)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c6)
Metrics::mape(mean(c6), conscious_value_for_c3_and_c6)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c6)
Metrics::smape(mean(c6), conscious_value_for_c3_and_c6)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c6)
fuzzyr.accuracy(mean(c6), conscious_value_for_c3_and_c6)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c3_and_c7)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c3_and_c7)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c3_and_c7)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c3_and_c8)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c3_and_c8)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c3_and_c8)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c3_and_c9)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c3_and_c9)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c3_and_c9)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c5)
Metrics::mape(mean(c5), conscious_value_for_c4_and_c5)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c5)
Metrics::smape(mean(c5), conscious_value_for_c4_and_c5)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c5)
fuzzyr.accuracy(mean(c5), conscious_value_for_c4_and_c5)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c6)
Metrics::mape(mean(c6), conscious_value_for_c4_and_c6)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c6)
Metrics::smape(mean(c6), conscious_value_for_c4_and_c6)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c6)
fuzzyr.accuracy(mean(c6), conscious_value_for_c4_and_c6)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c4_and_c7)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c4_and_c7)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c4_and_c7)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c4_and_c8)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c4_and_c8)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c4_and_c8)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c4_and_c9)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c4_and_c9)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c4_and_c9)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c6)
Metrics::mape(mean(c6), conscious_value_for_c5_and_c6)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c6)
Metrics::smape(mean(c6), conscious_value_for_c5_and_c6)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c6)
fuzzyr.accuracy(mean(c6), conscious_value_for_c5_and_c6)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c5_and_c7)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c5_and_c7)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c5_and_c7)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c5_and_c8)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c5_and_c8)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c5_and_c8)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c5_and_c9)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c5_and_c9)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c5_and_c9)

Metrics::mape(mean(c6), conscious_value_for_c6_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c6_and_c7)
Metrics::smape(mean(c6), conscious_value_for_c6_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c6_and_c7)
fuzzyr.accuracy(mean(c6), conscious_value_for_c6_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c6_and_c7)

Metrics::mape(mean(c6), conscious_value_for_c6_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c6_and_c8)
Metrics::smape(mean(c6), conscious_value_for_c6_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c6_and_c8)
fuzzyr.accuracy(mean(c6), conscious_value_for_c6_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c6_and_c8)

Metrics::mape(mean(c6), conscious_value_for_c6_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c6_and_c9)
Metrics::smape(mean(c6), conscious_value_for_c6_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c6_and_c9)
fuzzyr.accuracy(mean(c6), conscious_value_for_c6_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c6_and_c9)

Metrics::mape(mean(c7), conscious_value_for_c7_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c7_and_c8)
Metrics::smape(mean(c7), conscious_value_for_c7_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c7_and_c8)
fuzzyr.accuracy(mean(c7), conscious_value_for_c7_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c7_and_c8)

Metrics::mape(mean(c7), conscious_value_for_c7_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c7_and_c9)
Metrics::smape(mean(c7), conscious_value_for_c7_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c7_and_c9)
fuzzyr.accuracy(mean(c7), conscious_value_for_c7_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c7_and_c9)

Metrics::mape(mean(c8), conscious_value_for_c8_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c8_and_c9)
Metrics::smape(mean(c8), conscious_value_for_c8_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c8_and_c9)
fuzzyr.accuracy(mean(c8), conscious_value_for_c8_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c8_and_c9)

############################Humidity Validation Matrices########################################


Metrics::mape(mean(c1), conscious_value_for_c1_and_c2)
Metrics::mape(mean(c2), conscious_value_for_c1_and_c2)
Metrics::smape(mean(c1), conscious_value_for_c1_and_c2)
Metrics::smape(mean(c2), conscious_value_for_c1_and_c2)
fuzzyr.accuracy(mean(c1), conscious_value_for_c1_and_c2)
fuzzyr.accuracy(mean(c2), conscious_value_for_c1_and_c2)

Metrics::mape(mean(c2), conscious_value_for_c2_and_c3)
Metrics::mape(mean(c3), conscious_value_for_c2_and_c3)
Metrics::smape(mean(c2), conscious_value_for_c2_and_c3)
Metrics::smape(mean(c3), conscious_value_for_c2_and_c3)
fuzzyr.accuracy(mean(c2), conscious_value_for_c2_and_c3)
fuzzyr.accuracy(mean(c3), conscious_value_for_c2_and_c3)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c4)
Metrics::mape(mean(c4), conscious_value_for_c3_and_c4)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c4)
Metrics::smape(mean(c4), conscious_value_for_c3_and_c4)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c4)
fuzzyr.accuracy(mean(c4), conscious_value_for_c3_and_c4)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c5)
Metrics::mape(mean(c5), conscious_value_for_c3_and_c5)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c5)
Metrics::smape(mean(c5), conscious_value_for_c3_and_c5)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c5)
fuzzyr.accuracy(mean(c5), conscious_value_for_c3_and_c5)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c6)
Metrics::mape(mean(c6), conscious_value_for_c3_and_c6)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c6)
Metrics::smape(mean(c6), conscious_value_for_c3_and_c6)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c6)
fuzzyr.accuracy(mean(c6), conscious_value_for_c3_and_c6)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c3_and_c7)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c3_and_c7)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c3_and_c7)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c3_and_c8)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c3_and_c8)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c3_and_c8)

Metrics::mape(mean(c3), conscious_value_for_c3_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c3_and_c9)
Metrics::smape(mean(c3), conscious_value_for_c3_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c3_and_c9)
fuzzyr.accuracy(mean(c3), conscious_value_for_c3_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c3_and_c9)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c5)
Metrics::mape(mean(c5), conscious_value_for_c4_and_c5)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c5)
Metrics::smape(mean(c5), conscious_value_for_c4_and_c5)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c5)
fuzzyr.accuracy(mean(c5), conscious_value_for_c4_and_c5)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c6)
Metrics::mape(mean(c6), conscious_value_for_c4_and_c6)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c6)
Metrics::smape(mean(c6), conscious_value_for_c4_and_c6)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c6)
fuzzyr.accuracy(mean(c6), conscious_value_for_c4_and_c6)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c4_and_c7)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c4_and_c7)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c4_and_c7)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c4_and_c8)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c4_and_c8)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c4_and_c8)

Metrics::mape(mean(c4), conscious_value_for_c4_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c4_and_c9)
Metrics::smape(mean(c4), conscious_value_for_c4_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c4_and_c9)
fuzzyr.accuracy(mean(c4), conscious_value_for_c4_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c4_and_c9)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c6)
Metrics::mape(mean(c6), conscious_value_for_c5_and_c6)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c6)
Metrics::smape(mean(c6), conscious_value_for_c5_and_c6)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c6)
fuzzyr.accuracy(mean(c6), conscious_value_for_c5_and_c6)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c5_and_c7)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c5_and_c7)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c5_and_c7)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c5_and_c8)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c5_and_c8)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c5_and_c8)

Metrics::mape(mean(c5), conscious_value_for_c5_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c5_and_c9)
Metrics::smape(mean(c5), conscious_value_for_c5_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c5_and_c9)
fuzzyr.accuracy(mean(c5), conscious_value_for_c5_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c5_and_c9)

Metrics::mape(mean(c6), conscious_value_for_c6_and_c7)
Metrics::mape(mean(c7), conscious_value_for_c6_and_c7)
Metrics::smape(mean(c6), conscious_value_for_c6_and_c7)
Metrics::smape(mean(c7), conscious_value_for_c6_and_c7)
fuzzyr.accuracy(mean(c6), conscious_value_for_c6_and_c7)
fuzzyr.accuracy(mean(c7), conscious_value_for_c6_and_c7)

Metrics::mape(mean(c6), conscious_value_for_c6_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c6_and_c8)
Metrics::smape(mean(c6), conscious_value_for_c6_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c6_and_c8)
fuzzyr.accuracy(mean(c6), conscious_value_for_c6_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c6_and_c8)

Metrics::mape(mean(c6), conscious_value_for_c6_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c6_and_c9)
Metrics::smape(mean(c6), conscious_value_for_c6_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c6_and_c9)
fuzzyr.accuracy(mean(c6), conscious_value_for_c6_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c6_and_c9)

Metrics::mape(mean(c7), conscious_value_for_c7_and_c8)
Metrics::mape(mean(c8), conscious_value_for_c7_and_c8)
Metrics::smape(mean(c7), conscious_value_for_c7_and_c8)
Metrics::smape(mean(c8), conscious_value_for_c7_and_c8)
fuzzyr.accuracy(mean(c7), conscious_value_for_c7_and_c8)
fuzzyr.accuracy(mean(c8), conscious_value_for_c7_and_c8)

Metrics::mape(mean(c7), conscious_value_for_c7_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c7_and_c9)
Metrics::smape(mean(c7), conscious_value_for_c7_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c7_and_c9)
fuzzyr.accuracy(mean(c7), conscious_value_for_c7_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c7_and_c9)

Metrics::mape(mean(c8), conscious_value_for_c8_and_c9)
Metrics::mape(mean(c9), conscious_value_for_c8_and_c9)
Metrics::smape(mean(c8), conscious_value_for_c8_and_c9)
Metrics::smape(mean(c9), conscious_value_for_c8_and_c9)
fuzzyr.accuracy(mean(c8), conscious_value_for_c8_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c8_and_c9)

############################Validation Matrices########################################
mape(mean(c8), conscious_value_for_c8_and_c9)
mape(mean(c9), conscious_value_for_c8_and_c9)
smape(mean(c8), conscious_value_for_c8_and_c9)
smape(mean(c9), conscious_value_for_c8_and_c9)
fuzzyr.accuracy(mean(c8), conscious_value_for_c8_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c8_and_c9)

mape(mean(c8), conscious_value_for_c8_and_c9)
mape(mean(c9), conscious_value_for_c8_and_c9)
smape(mean(c8), conscious_value_for_c8_and_c9)
smape(mean(c9), conscious_value_for_c8_and_c9)
fuzzyr.accuracy(mean(c8), conscious_value_for_c8_and_c9)
fuzzyr.accuracy(mean(c9), conscious_value_for_c8_and_c9)



