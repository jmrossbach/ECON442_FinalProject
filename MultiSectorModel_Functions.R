#################################################################
##This File Includes Functions that Will be Used to Solve Model##
#################################################################
##Create a function for turning vectors into Arrays
vec2array <- function(X,dim,size_vec=NULL) {
  ##Optionally, can provide size of each dimension, otherwise assume MxJxJ Array
  if(!is.null(size_vec)) {
    r = size_vec[1]; c=size_vec[2]; m=size_vec[3];
  } else {
    r = M; c=J; m=J;
  }
  
  if(dim==1) { ##If dim is 1, then dimension is sector
    result <- array(rep(matrix(rep(X,c), ncol=c, byrow=FALSE),m),dim=c(r,c,m))
    
  } else if(dim==2) { ##dim==2 means dimension is destination
    result <- array(rep(matrix(rep(X,r), nrow=r, byrow=TRUE),m),dim=c(r,c,m))
    
  } else if(dim==3) { ##dim==3 means dimension is source
    
    temp_vector <- vector()
    
    for(j in 1:J) {
      temp_vector <- c(temp_vector,rep(X[j],r*c))
    }
    result <- array(temp_vector,dim=c(r,c,m))
  }
  
  
  return(result)
}


##Create a function for turning vectors into Arrays
mat2array <- function(X,missing_dim,size_vec=NULL) {
  ##Optionally, can provide size of each dimension, otherwise assume MxJxJ Array
  if(!is.null(size_vec)) {
    r = size_vec[1]; c=size_vec[2]; m=size_vec[3];
  } else {
    r = M; c=J; m=J;
  }
  
  if(missing_dim==1) { ##If dim is 1, then dimension is sector
    #result <- array(rep(X,r),dim=c(r,c,m))
    temp <- array(rep(X,r),dim=c(c,r,m))
    result <- array(,dim=c(r,c,m)) ##Purposefully empty
    for(i in 1:m) {
      result[,,i] <- t(temp[,,i])
    }
  } else if(missing_dim==2) { ##dim==2 means dimension is destination
    result <- array(rep(X,c),dim=c(r,c,m))
  } else if(missing_dim==3) { ##dim==3 means dimension is source
    result <- array(rep(X,m),dim=c(r,c,m))
  }
  
  return(result)
  
}

##########################################################################
##Functions for Reshaping Endogenous Variables into Vector for rootSolve##
##########################################################################
vars2vec <- function(l,y,p,P_index,Y_index,T,Income) {
  vec_list <- vector()
  vec_list <- c(vec_list,as.vector(l))
  vec_list <- c(vec_list,as.vector(y))
  vec_list <- c(vec_list,as.vector(p))
  vec_list <- c(vec_list,as.vector(P_index))
  vec_list <- c(vec_list,as.vector(Y_index))
  vec_list <- c(vec_list,as.vector(T))
  vec_list <- c(vec_list,as.vector(Income))
  
  return(vec_list)
}

vec2vars <- function(vec_list) {
  index = 1
    #w <- vec_list[index:(J-1)]
    #index = index+(J-1)
  ##MxJxJ Arrays

  l <- array(vec_list[index:(index+M*J*J-1)],dim<-c(M,J,J))
  index = index+(M*J*J)
  y <- array(vec_list[index:(index+M*J*J-1)],dim<-c(M,J,J))
  index = index+(M*J*J)
  p <- array(vec_list[index:(index+M*J*J-1)],dim<-c(M,J,J))
  index = index+(M*J*J)
  
  ##MxJ Matrices
  P_index <- matrix(vec_list[index:(index+M*J-1)],nrow=M,ncol=J)
  index = index+(M*J)
  Y_index <- matrix(vec_list[index:(index+M*J-1)],nrow=M,ncol=J)
  index = index+(M*J)
  
  ##Jx1 vectors
  T <- vec_list[index:(index+J-1)]
  index = index+(J)
  Income <- vec_list[index:(index+J-1)]
  index = index+(J)

  #varlist <- list("w"=w,"l"=l,"y"=y,"p"=p)
  varlist <- list("l"=l,"y"=y,"p"=p,"P_index"=P_index,"Y_index"=Y_index,"T"=T,"Income"=Income)
  return(varlist)
}


#############################################################
##WRITE FUNCTION FOR EQUILIBRIUM
#############################################################
equilibrium_eqns <- function(w_guess,other_vars) {
	W <- c(1,w_guess) ##Walras law, so normalize wage of first country to 1

	########################################################
	##Endogeneous Parameters of the Model
	########################################################
  vars<- vec2vars(other_vars)
  l <- vars$l
  y <- vars$y
  p <- vars$p
  P_index <- vars$P_index
  Y_index <- vars$Y_index
  T <- vars$T
  Income <- vars$Income
  
	########################################################
	##Equilibrium Conditions
	#########################################################
  
  ##We really only need to solve for wages, so only have labor clearing as our test
  ##For the rest we don't need rootSolve, we can simply iterate and they will converge.  Much faster than newton's method on a giant matrix.
  max_diff = 1; ##Keep track of how far the variables are from their previous iteration value
  while(max_diff>1e-6) {
    vec_old <- vars2vec(l,y,p,P_index,Y_index,T,Income) ##Store old values, for doing diff

    p <- (tau/z)/vec2array(W,3) ##Update p with Firm's Problem
    P_index <- rowSums(mu^vec2array(sigma,1)*p^vec2array(1-sigma,1),dims=2)^(1/(1-sigma)) ##P index is definition of price index
    
    y <- (mu/p)^vec2array(sigma,1)*mat2array(Y_index,3)*(mat2array(P_index,3)^(vec2array(sigma,1))) ##Update y with solution to CP
    l <- d*y/z ##Update labor with production function
    
    T<-apply((tau-1)*p*y,2,sum) ##Tariff Revenue.  Need to recompute with new trade flows and tariffs (note tau, not tau_before)
    Income <- W*L + T + D ##Income is definition of income.  Wages change, Tariffs change.  L and D are exogenous.
    
    for(j in 1:J) {
      Y_index[,j] = (theta[,j]/P_index[,j])*Income[j] ##Y index is solution to CP also
    }
    
    vec_new <- vars2vec(l,y,p,P_index,Y_index,T,Income) ##Store updated values, for doing diff
  	max_diff=max(abs(vec_new-vec_old)) ##see how far each one is from the old version
  }

  ##Labor Market Clearing Conditions
	L_Diff <- L - apply(l,3,sum); ##Sum l values over source (3rd dim)
  L_Diff <- L_Diff[2:J] ##Walras' Law, ignore first country
	
  ##   
  return(L_Diff)  ##Vector showing how far each labor market equation is from binding

}

##########################################################################
##Functions for Assigning Global Variables##
##########################################################################
assign_equilibrium <- function(W,suffix=NULL) {
  max_diff = 1; ##Keep track of how far the variables are from their previous iteration value
  while(max_diff>1e-6) {
    vec_old <- vars2vec(l,y,p,P_index,Y_index,T,Income) ##Store old values, for doing diff
    
    p <- (tau*d/z)/vec2array(W,3) ##Update p with Firm's Problem
    P_index <- rowSums(mu^vec2array(sigma,1)*p^vec2array(1-sigma,1),dims=2)^(1/(1-sigma)) ##P index is definition of price index
    
    y <- (mu/p)^vec2array(sigma,1)*mat2array(Y_index,3)*(mat2array(P_index,3)^(vec2array(sigma,1))) ##Update y with solution to CP
    l <- d*y/z ##Update labor with production function
    
    T<-apply((tau-1)*value_before,2,sum) ##Tariff Revenue.  Need to recompute with new trade flows and tariffs (note tau, not tau_before)
    Income <- W*L + T + D ##Income is definition of income.  Wages change, Tariffs change.  L and D are exogenous.
    
    for(j in 1:J) {
      Y_index[,j] = (theta[,j]/P_index[,j])*Income[j] ##Y index is solution to CP also
    }
    
    vec_new <- vars2vec(l,y,p,P_index,Y_index,T,Income) ##Store updated values, for doing diff
    max_diff=max(abs(vec_new-vec_old)) ##see how far each one is from the old version
  }
  
  ##Assign the variables, with any suffix
  assign(paste("p",suffix,sep = ""), p, envir = .GlobalEnv)
  assign(paste("y",suffix,sep = ""), y, envir = .GlobalEnv)
  assign(paste("c",suffix,sep = ""), y, envir = .GlobalEnv)
  assign(paste("l",suffix,sep = ""), l, envir = .GlobalEnv)
  assign(paste("T",suffix,sep = ""), T, envir = .GlobalEnv)
  assign(paste("Income",suffix,sep = ""), Income, envir = .GlobalEnv)
  assign(paste("P_index",suffix,sep = ""), P_index, envir = .GlobalEnv)
  assign(paste("Y_index",suffix,sep = ""), P_index, envir = .GlobalEnv)
}



