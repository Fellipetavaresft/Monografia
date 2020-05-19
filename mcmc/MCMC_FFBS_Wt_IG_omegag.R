#####################################################################################################################
#	 	FFBS para MLD Hierárquico Normal Multivariado em 4 estágios									
#					
#####################################################################################################################
				
#####################################################################################################################
#			Início do Programa	
#####################################################################################################################

rm(list=ls(all=TRUE))

library(mvtnorm)
library(magic)
library(MCMCpack)

##############################
################################################
#		MCMC + Forward Filtering#
##############################################################


for(it in 1:totit)
{
	# Updating: Atualização de t = 1.
	t=1
		#Priori em t=1
		at[t] = Gt[t]*m0
		Rt[t] = ((Gt[t])^2)*C0 + Wt[t]

		#Distribuição Preditiva
		ft[t,] = F4t[,,t]*at[t]
		Qt[t,,] = (F4t[,,t]*Rt[t])%*%t(F4t[,,t]) + V4t[,,t]

		#Posteriori em 1
		et[t,] = (theta3t[t,] - ft[t,])
		At[t,,] = (Rt[t]*t(F4t[,,t]))%*%solve(Qt[t,,])						   
   		mt[t] = at[t] + At[t,,]%*%et[t,]
    		Ct[t] = Rt[t] - (t(At[t,,])%*%Qt[t,,]%*%At[t,,])
		
	# Updating: Atualização de t = 2, ... ,(tstar-1)
	for(t in 2:T)
	{
		#Priori em t=1
		at[t] = Gt[t]*mt[t-1]
		Rt[t] = ((Gt[t])^2)*Ct[t-1] + Wt[t]
	
		#Distribuição Preditiva
		ft[t,] = F4t[,,t]*at[t]
		Qt[t,,] = (F4t[,,t]*Rt[t])%*%t(F4t[,,t]) + V4t[,,t]

		#Posteriori em 1
		et[t,] = (theta3t[t,] - ft[t,])
		At[t,,] = Rt[t]*t(F4t[,,t])%*%solve(Qt[t,,])						   
   		mt[t] = at[t] + At[t,,]%*%et[t,]
   		Ct[t] = Rt[t] - (t(At[t,,])%*%Qt[t,,]%*%At[t,,])
	}

	#Smoothing T,...,1
	#em T
		ms[T] = mt[T]
		Cs[T] = Ct[T]
		theta4t[T] = rnorm(1,ms[T],sqrt(Cs[T]))

	#em t=T-1,...,tstar+1
	for(t in (T-1):1)
	{
		ms[t] = mt[t] + (Ct[t]%*%t(Gt[t+1])%*%solve(Rt[t+1])%*%(theta4t[t+1]-Gt[t+1]%*%mt[t]))
		Cs[t] = Ct[t] - (Ct[t]%*%t(Gt[t+1])%*%t(solve(Rt[t+1]))%*%(Gt[t+1])%*%t(Ct[t]))
		theta4t[t,] = rnorm(1,ms[t],sqrt(Cs[t]))
	}

	##Sorteando Wt
	auxWt = 0
	auxWt = ((theta4t[1] - theta4_0 )^2)
	for(t in 2:T)
	{
		auxWt = auxWt + ((theta4t[t] - theta4t[t-1])^2)
	}
	Wt[1] = (1/rgamma(1,alphaW +(T/2),betaW + ((1/2)*auxWt)))
	
	for(t in 2:T)
	{
		Wt[t] = Wt[1]
	}

	##Sorteando Tau2
	somatau = 0
	for(t in 1:T)
	{
		cont=1
		for (g in 1:G)
		{
			for (i in 1:Ng[g])
			{	
				for(j in 1:J)
				{
					somatau = somatau + ((Yt[t,cont]-theta1t[t,cont])^2)
					cont = cont + 1
				}
			}
		}
	}
	tau2 = (1/rgamma(1,atau+(I*J*T/2),btau+((1/2)*somatau)))
	
	##Atualizando as cadeias

	if ((it>=burnin)&&((it%%thin)==0))	   	
	{
		v.tau2[conta] = tau2
		v.Wt[conta] = Wt[1]
				for (t in 1:T)
		{
			v.theta1t[conta ,t,] = theta1t[t,]
			v.theta2t[conta ,t,] = theta2t[t,]
			v.theta3t[conta ,t,] = theta3t[t,]
			v.theta4t[conta ,t] = theta4t[t]
		}
		conta = conta + 1
	}
	print(it)
}
