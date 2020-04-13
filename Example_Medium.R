#############################################################
##INSTALL AND LOAD EXTERNAL R PACKAGES
#############################################################
##Install packages "rootSolve" and "ggplot2", .,
#install.packages("rootSolve")  ##Install rootSolve if you don't already have it installed.  Comment this line out if you already have the package.
#install.packages("ggplot2") ##Install ggplot2 if you don't already have it installed.  Comment this line out if you already have the package.
#install.packages("dplyr") ##Install dplyr if you don't already have it installed.  Comment this line out if you already have the package.
#install.packages("openxlsx") ##Install openxlsx if you don't already have it installed.  Comment this line out if you already have the package.
rm(list = ls()) ##Clear the workspace
library("rootSolve")  ##Load the library, since we will be using the function multiroot()
library("ggplot2") ##For Plotting data
library("dplyr") ##For Plotting data
library("openxlsx") ##For exporting data to excel
source('MultiSectorModel_Functions.R') ##Load the functions we will use to solve the model

######################################
##STEP 1 LOAD DATA.  Store in Arrays##
######################################
## Before running code go to the menu in RSTudio and select:
## Session -> Set Working Directory -> To Source File Location
WIOD_data <- read.csv("Medium_Size_Data.csv",header=TRUE)  ##Read in Data.  Change name of csv file if different file.

J <- length(unique(WIOD_data$source))  ##Number of countries
M <- length(unique(WIOD_data$sector))  ##Number of sectors

##Store variables in arrays.  Row (first dimension) is sector, Column (second dimension) is destination country, Matrix (third dimension) is source country
##This means we need to sort first by source, then dest, then sector
WIOD_data <- arrange(WIOD_data,source,destination,sector)

##Names for our Arrays, if we want to attach them.  I don't want to right now, so comment out.
##To include names, do ", dimnames = list(row.names,col.names,mat.names)" at the end of the array() function
row.names <- unique(WIOD_data$sector)
col.names <- unique(WIOD_data$destination)
mat.names <- unique(WIOD_data$source)
name_list <- list(row.names,col.names,mat.names)
index = paste(WIOD_data$source, WIOD_data$sector,WIOD_data$destination,sep="-")
WIOD_data$index <- index

##BASE YEAR VALUES  ##Replace 2000 with your BASE YEAR (up to 2015)
value_before <- array(WIOD_data$value2000,dim=c(M,J,J), dimnames=name_list)  ##Trade Flows in Base Year
tau_before <- array(WIOD_data$simpleaverage2000,dim=c(M,J,J), dimnames=name_list)  ##Tariffs in Base Year
Z_before <- array(WIOD_data$Z_2000,dim=c(M,J,J), dimnames=name_list)  ##Sector Productivity in Base Year
iceberg_before <- array(rep(1,M*J*J),dim=c(M,J,J), dimnames=name_list) ##Iceberg costs in Base Year, can assume 1 if don't have data

##AFTER VALUES (Don't need anythign for value_after, do for the rest).  You can take from Data or just change yourself
value_after <- array(WIOD_data$value2006,dim=c(M,J,J), dimnames=name_list)  ##Trade Flows in 2006 (Includes self-trade)
tau_after <- array(WIOD_data$simpleaverage2006,dim=c(M,J,J), dimnames=name_list)  ##Tariffs in 2006
iceberg_after <- array(rep(1,M*J*J),dim=c(M,J,J), dimnames=name_list) ##Iceberg 
Z_after <- array(WIOD_data$Z_2006,dim=c(M,J,J), dimnames=name_list)  ##Change in Sector Productivity from 2000 to 2006

##Suppose you want to do a hypothetical counterfactual, and you want to change tariffs for a given importer, exporter, sector pair.  Then can do something like follows
#tau_after <- tau_before ##Set tariffs to base value
#tau_after[WIOD_data$sector == "Agriculture" & WIOD_data$source == "CHN"] <- 1.15*tau_before[WIOD_data$sector == "Agriculture" & WIOD_data$source == "CHN"]  ##Increase them only for exports in Agriculture from Chian


###############################
##STEP 2 Calibrate Model.    ##
###############################
##First, need Elasticities from outside of the model.
sigma <- rep(8,M) ##Could also just have it equal to 8, but doing this so we could have different elasticities for different sectors later
#sigma <- array(rep(8,M*J*J),dim=c(M,J,J))  

#Results <- function(value_before,tau_before,sigma) {
W <- rep(1,J)  ##All wages equal to 1
z <- array(rep(1,M*J*J),dim=c(M,J,J))  ##All productivities equal to 1

Expenditures <- value_before * tau_before ##Trade values do not include tariffs, so add them
##For tariffs and Trade Balance, use Trade Values instead of Expenditures  
T<-apply((tau_before-1)*value_before,2,sum) ##Tariff Revenue.  Sum for each destination (2nd dim) 
D<-apply(value_before,2,sum) - apply(value_before,3,sum) ##Trade Deficit.  Compare total trade as dest (2nd dim) vs total trade as source (3rd dim)
##For GDP/Income, use Expenditures
Income <- apply(Expenditures,2,sum)  #Income.  Sum all expenditure for destination (includes self trade)
L <- (Income-T-D)/W  ##Labor supply is Total income, minus trade balance and tariff income, divided by wage

##find l,p,y
l <- value_before/vec2array(W,3) ###Labor is zero profit condition for firm.  Use the wage of source country, dim=3
p <- (tau_before*iceberg_before/z)/vec2array(W,3) ##p is from the sol to firm's prob.  Use wage of source country, dim=3
y <- (Expenditures/p) ##y is simply from definition of expenditures = p*y [note c=y]

mu <- (p*y^vec2array(1/sigma,1)) / mat2array(rowSums(p*y^vec2array(1/sigma,1),dims=2),3)  ##Mu is from calibration formula

P_index <- rowSums(mu^vec2array(sigma,1)*p^vec2array(1-sigma,1),dims=2)^(1/(1-sigma)) ##P index is definition of price index
Y_index <- rowSums(Expenditures,dims=2)/P_index

theta <- matrix(,nrow=M,ncol=J)  ##Purposefully empty
for(j in 1:J) {
  theta[,j] = (P_index[,j]*Y_index[,j])/Income[j]
}

########################################################
##Step 2a: Save pre-Counterfactual Values for Later Comparison
########################################################
tau = tau_before
d = iceberg_before
assign_equilibrium(W,suffix="_old") ##Assign variables for the Old Equilibrium, each variable has suffix _New

#############################################################
##Step 3: Put in CounterFactuals
#############################################################
other_vars <- vars2vec(l,y,p,P_index,Y_index,T,Income); ##other_vars is guess for all non-wage endogenous variables
w_guess <- W[2:J] ##w_guess is our guess for 

##change in iceberg costs (but make sure self trade stays 1)
delta_iceberg <- iceberg_after/iceberg_before ##Change in Iceberg Costs from 2000 to 2006, .9 corresponds to 10 percent decrease; 1 means no change
for(i in 1:J) {delta_iceberg[,i,i] <- 1} #Iceberg costs don't change for domestic trade, since always 1

##Change Tariffs and Sectoral Productivities
tau = tau_after
d = iceberg_before * delta_iceberg
z = z*mat2array(Z_after/Z_before,missing_dim = 2)

########################################################
##Step 4: Run Nonlinear Solver to Find New Equilibrium##
########################################################
solution <- multiroot(f=equilibrium_eqns,start=w_guess,parms=other_vars)  ##type "help(multiroot)" without quotes to learn more about this function
solution ##See whether it worked.  Look at tolerance for how far you are from zero. If far from zero then didn't converge and need a better initial guess.

W_new <- c(1,solution$root) ##multiroot solves for the new wage vector
assign_equilibrium(W_new,suffix="_new") ##Assign variables for the New Equilibrium, each variable has suffix _New

#############################
###Step 5: Explore Results###
#############################
###Compute some percent changes.  Look at trade divided by GDP for everything.
Trade_After_Model = p_new*y_new/Income_new
Trade_After_Data = value_after/Income_new
TradeChanges_Model = 100*((p_new*y_new/Income_new)-(p_old*y_old/Income_old))/(p_old*y_old/Income_old)
TradeChanges_Data = 100*(value_after/Income_new-value_before/Income_old)/(value_before/Income_old)  ##Only have this line if you have changes in data before and after that you want

sprintf("Correlation between Model and Data is %s",cor(as.vector(TradeChanges_Model),as.vector(TradeChanges_Data),use="complete.obs")) #complete obs means exclude nans

##Add results to our table
Results <- WIOD_data
Results$TradeChanges_Model <- as.vector(TradeChanges_Model) 
Results$TradeChanges_Data <- as.vector(TradeChanges_Data)
Results$Trade_After_Model <- as.vector(Trade_After_Model) 
Results$Trade_After_Data <- as.vector(Trade_After_Data)

##Export them to a excel file so you can make graphs in Excel
write.xlsx(Results, file="Model_Results.xlsx")
##If above command doesn't work, then just write to a csv file and open in excel yourself
#write.csv(Results, file="Model_Results.csv")

##############################################
##Step 6: Make Graphs in RStudio (Optional)
##############################################
##Plot Changes in Data versus Model
plot(Results$TradeChanges_Data[Results$TradeChanges_Data<300 & Results$TradeChanges_Model < 300],Results$TradeChanges_Model[Results$TradeChanges_Data<300 & Results$TradeChanges_Model < 300])
##Use ggplot to add a regression line (requires loading ggplot2 library)
ggplot(Results,aes(x=TradeChanges_Data, y=TradeChanges_Model)) + geom_point(shape=1) + geom_smooth(method=lm) #ggplot2
ggplot(Results,aes(x=TradeChanges_Data, y=TradeChanges_Model)) + scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2')  + geom_point(shape=1) + geom_smooth(method=lm) #Log Scale

##Use ggplot to add a regression line (requires loading ggplot2 library)
plot(Results$Trade_After_Data,Results$Trade_After_Model, log="xy") ##Use a log scale
ggplot(Results,aes(x=Trade_After_Data, y=Trade_After_Model,label=index)) + scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2') + geom_point(shape=1) + geom_text(aes(label=index))+ geom_smooth(method=lm) #Add labels


###LOOK AT CHN Exports
CHN_USA_Results <- subset(Results,source=="CHN" & destination=="USA") ##Focus on Only China's Exports to USA
##Graph for Exports in Model
ggplot(CHN_USA_Results, aes(x=sector,y=TradeChanges_Model,fill=destination))  + theme(axis.text.x=element_text(angle=90,hjust=1)) +geom_bar(position="dodge",stat="identity")
ggplot(CHN_USA_Results, aes(x=sector,y=TradeChanges_Data,fill=destination))  + theme(axis.text.x=element_text(angle=90,hjust=1)) +geom_bar(position="dodge",stat="identity")
##Regression of changes
ggplot(CHN_USA_Results,aes(x=TradeChanges_Data, y=TradeChanges_Model)) + geom_point(shape=1) + geom_smooth(method=lm) #ggplot2