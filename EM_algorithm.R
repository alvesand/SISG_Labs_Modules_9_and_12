#SISG Module 12 - Lab3: Programming the EM algorithm for Variance Component Estimation
#Packages required for this class 
#install.packages("Rfast")

#Lets create our own function for variance component estimation using the EM algorithm (Just a few lines of code)
emreml610 = function(maxrounds, guess, criter, Ainv, X, Z, y){
  library(MASS) 
  library(Rfast) #Will be used for faster matrix multiplications
  n = length(y) #number of non missing observations
  p=ncol(X) #number of effects in the X matrix
  q=ncol(Z) #number of animals in the pedigree file (allows q > n)
  
  old_var_u = guess[1]
  old_var_e = guess[2]
  lambda = as.numeric(old_var_e/old_var_u)
  rankX = qr(X)$rank #Number of independent columns in X
  
  #Build all MME unchanged blocks 
  B11 = crossprod(X) #same as t(X)%*%X, but a little bit faster
  B12 = crossprod(X,Z)
  B21 = crossprod(Z,X)
  ZZ = crossprod(Z,Z)
  LHS1 = cbind(B11,B12) #The upper block of the Left-hand side
  RHS = rbind(crossprod(X,y),crossprod(Z,y))
  
  ini = Sys.time()
  #Store the user starting guesses
  guesses = matrix(0,maxrounds,2)
  guesses[1,] = guess
  cat('Start EM interations!') #A little text art for the verbose to look fancy =)
  cat('\n')
  cat('
     /)  (\
.-._((.~~.))_.-.
 `-.   @@   .-.
   / .o--o. \
  ( ( .__. ) )
   
   ANSCI610'
  ) 
  cat('\n')
  
  for (i in 2:maxrounds){
    Gi = Ainv*as.numeric(lambda) 
    B22 = ZZ+Gi
    
    LHS = rbind(LHS1,cbind(B21,B22)) #Complete the LHS
    LHSi = chol2inv(chol(LHS)) #Like solve() but faster!
    
    C22 = LHSi[(p+1):nrow(LHS),(p+1):nrow(LHS)]
    
    #mat.mult(A,B) is equivalent to A%*%B but a little bit faster (Require Rfast)
    sol = mat.mult(LHSi,RHS) #Get the Mixed Model Solutions
    beta = as.matrix(sol[1:p,1]) #Separate the BLUE
    u = as.matrix(sol[(p+1):(q+p),1]) #Separate the BLUP
    
    yhat = mat.mult (X,beta)+mat.mult(Z[,(q-n+1):q],as.matrix(u[(q-n+1):q])) 
    e = y - yhat #Get the residuals
    
    var_e = as.numeric(crossprod(e,y)/(n-rankX))
    var_u = (mat.mult(crossprod(u,Ainv),u)+sum(Ainv*C22)*var_e)/q
    guesses[i,] = c(var_u,var_e)
    #Print the current iteration
    cat('\n')
    cat("Round: ", i-1, '\n')
    cat("G", guesses[i,1], '\n')
    cat("R", guesses[i,2], '\n')
    lambda = as.numeric(var_e/var_u) #Update lambda
    #Assess the algorithm stopping criteria
    if (max(abs(guesses[i,]-guesses[i-1,]))<criter){
      break
    }
  }
  
  #Compute final results
  PEV = c(diag(LHSi*var_e))
  blue = data.frame("beta" = beta, s.e = sqrt(PEV[1:p]))
  blup = data.frame("ebv" = u, s.e = sqrt(PEV[(p+1):(q+p)]), acc = sqrt(1-PEV[(p+1):(q+p)]/c(var_u)))
  fim = Sys.time()
  
  cat('\n')
  cat("Time elapsed: ", fim-ini, '\n')
  return(list(G = guesses[i,1],R = guesses[i,2], blue= blue, blup = blup))
}

