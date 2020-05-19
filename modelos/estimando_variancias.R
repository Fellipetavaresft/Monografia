rm(list=ls(all=TRUE))

library("dlm")
library("forecast")
library("tseries")
#library("ggplot2")
#library("gridExtra")
#library("grDevices")
#library("grid")
#library("tidyverse")

dados<-read.table("IPCA.txt",header=T,sep=",")
attach(dados)

#dados corrigidos para o periodo correto 1995 - 2015
ipca<-c(IPCA[181:432])
ts.ipca=ts(ipca,frequency=12,start=c(1995)) 

#############################################################################################
#					DLM utilizado no Trabalho
#############################################################################################

m = rep(0,12)
c = 1e-6 * diag(12) 
F = matrix(c(1,1,0,0,0,0,0,0,0,0,0,0), nr=1)

G = bdiag(1,matrix(c(
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
1,0,0,0,0,0,0,0,0,0,0,
0,1,0,0,0,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,0,0,0,
0,0,0,1,0,0,0,0,0,0,0,
0,0,0,0,1,0,0,0,0,0,0,
0,0,0,0,0,1,0,0,0,0,0,
0,0,0,0,0,0,1,0,0,0,0,
0,0,0,0,0,0,0,1,0,0,0,
0,0,0,0,0,0,0,0,1,0,0,
0,0,0,0,0,0,0,0,0,1,0),nr=11,byrow=TRUE))
V = 0.003
#W = diag(c(20000, 400, 0, 0,0,0,0,0,0,0,0,0,0), nr=12)
#W = diag(c(2000, 4000, 100, 100,100,100,0,0,0,0,0,0,0), nr=12)

W = diag(c(1,1, 0, 0,0,0,0,0,0,0,0,0,0), nr=12)


model.1 = dlm(m0 = m, C0 = c, FF = F , GG = G, V = V, W = W )
output.1 = dlmFilter(ts.ipca,model.1)
Suaviz.1 = dlmSmooth(output.1)


plot(ts.ipca, type="p", ylab="Expenditure")
lines(output.1$m[,1],type="l", col = "red")

#############################################################################################
#					Estimando variâncias - Gibbs
#############################################################################################

set.seed(5672)
time_ini = Sys.time()
MCMC <- 20000 #20000
gibbsOut.2 <- dlmGibbsDIG(ts.ipca, mod = model.1, a.y = 1,b.y = 1000, a.theta = 1, b.theta = 1000, n.sample = MCMC,thin = 0, save.states = FALSE)
#gibbsOut.2 <- dlmGibbsDIG(ts.ipca, mod = model.1,a.y = 1,b.y = 4, a.theta =1, b.theta = .00009, n.sample = MCMC,thin = 0, save.states = FALSE)
###gibbsOut.2 <- dlmGibbsDIG(ts.ipca, mod = model.1,shape.y = .01,rate.y = .01, shape.theta = 10, rate.theta = .01, n.sample = MCMC,thin = 0, save.states = FALSE)
#gibbsOut.2 <- dlmGibbsDIG(ts.ipca, mod = model.1,shape.y = .01,rate.y = .01, shape.theta = 2, rate.theta =1, n.sample = MCMC,thin = 0, save.states = FALSE)

#gibbsOut.2 <- dlmGibbsDIG(ts.ipca, mod = model.1,shape.y = .01,rate.y = .01, shape.theta = 0.1, rate.theta = .01, n.sample = MCMC,thin = 0, save.states = FALSE)
time_final = Sys.time()

# testar #
#gibbsOut.2 <- dlmGibbsDIG(ts.ipca, mod = model.1,shape.y = .1,rate.y = .1, shape.theta = 0.1, rate.theta = .1, n.sample = MCMC,thin = 0, save.states = FALSE)


#a.y=.01 ,  b.y=.01, a.theta=2, b.theta=1
#burn = 2001
#use <- MCMC - burn
#from <- 0.05 * use
#thin = 10
#dV[seq(burn:MCMC:thin)]

burn = 1000
thin = 10
length(seq(burn,MCMC,thin))
#cadeia = seq(burn,MCMC,thin)
cadeia = seq(burn,MCMC,thin)


plot(gibbsOut.2$dW[cadeia,1],type="l")

hist(gibbsOut.2$dW[,5])
gibbsOut.2$dW
plot(ergMean(gibbsOut.2$dV[-(1:burn)], from), type="l", xaxt="n",xlab="", ylab="")
at <- pretty(c(0,use),n=3); at <- at[at>=from]
axis(1, at=at-from, labels=format(at))

mcmcMeans(cbind(gibbsOut.2$dV[-(1:burn)], gibbsOut.2$dW[-(1:burn),]))

#matriz_w <- mcmcMeans(gibbsOut.2$dW[2001,])
#matriz_v <- mcmcMeans(gibbsOut.2$dV[-(1:burn)]) 
var(ts.ipca)
matriz_w <- mcmcMeans(gibbsOut.2$dW[seq(burn,MCMC,thin),])
matriz_v <- mcmcMeans(gibbsOut.2$dV[seq(burn,MCMC,thin)]) 
matriz_w
matriz_v
#############################################################################################
#					Modelo com novas variâncias
#############################################################################################

W = diag(c(matriz_w[1,]), nr=12)
V = matriz_v[1]

model.2 = dlm(m0 = m, C0 = c, FF = F , GG = G, V = V, W = W )
output.2 = dlmFilter(ts.ipca,model.2)
Suaviz.2 = dlmSmooth(output.2)

###############################################################################

set.seed(5672)
nfuturo = 1
fore <- dlmForecast(output.2,nAhead = nfuturo)

attach(fore)
f = fore$f


k<-253:264
x<-seq(0,length(ipca),1)
attach(fore)

hwidth = qnorm(0.05/2, lower = FALSE) * sqrt(unlist(fore$Q))
pl <- fore$a[,1] + qnorm(0.05, lower = FALSE) * sqrt(unlist(fore$Q))

pu <- fore$a[,1] + qnorm(0.95, lower = FALSE) * sqrt(unlist(fore$Q))

windows(width = 5.9, height = 4, rescale = "R")
par(mar=c(4,4,2,0))#bottom,left,top,rigth#fore2 <- cbind(f, as.vector(f) + hwidth %o% c(-1, 1))
#x11(5.90,2.43)
plot(ipca, xlab="", ylab="",type="p",main="",pch=3,xlim=c(0,276), ylim=c(),lwd=1 ,cex.lab=1,cex=1,cex.axis=1,axes=FALSE)
title(ylab=list(expression("IPCA Variação Mensal")),cex.lab=1)
title(xlab=list(expression("Período")),cex.lab=1)
axis(2,cex.axis=1)
axis(1,seq(1,276,12),(1995:2017),cex.axis=1)
#lines(x,output.2$m[,1],lty=1,lwd=2,col="green") #ajuste hw
lines(x,Suaviz.2$s[,1],lty=1,lwd=2) #ajuste hw
lines(k,pl,lty = 2,lwd=2)
lines(k,pu,lty = 2,lwd=2)
abline(v=252,lty=2)
lines(k,IPCA[433:444],type="p",cex.lab=1.4,pch=19)
lines(k,f,lty=1,lwd=2,col="blue") # previsao
#points(fore$newObs[[2]], pch = 16, col = "6")


#############################################################################################
#					Comparações
#############################################################################################

N = length(output.2$m[,1])
n = length(ts.ipca)
#MSE

MSE = 1/length(ts.ipca)*(sum((ts.ipca - output.2$m[2:N,1])^2));MSE

MAPE = ((sum(abs((ts.ipca - output.2$m[2:N,1])/(ts.ipca))))/length(ts.ipca))*100;MAPE

SMAPE = ((sum(abs((ts.ipca - output.2$m[2:N,1]))/((ts.ipca + output.2$m[2:N,1])/2)))/length(ts.ipca))*100;SMAPE


U =sqrt( (sum((ts.ipca - output.2$m[2:N,1])^2))/(sum((ts.ipca[2:n] - ts.ipca[1:n-1])^2)) );U

MAD = (sum(abs(ts.ipca - output.2$m[2:N,1])))/(length(ts.ipca));MAD


mle2 <- dlmMLE(y = ts.ipca, parm = rep(0,12), build = model2)
aic2 <- -2 * mle2$value + 2 * 12

# Data frame com resultados das métricas
resultado = data.frame(MSE,MAPE,U)
# Gerar uma tabela que será visualizada no documento. Mais detalhes em http://haozhu233.github.io/kableExtra/
knitr::kable(resultado, format = "pandoc", digits = c(2,2,2), align = 'c')


##############################################################################

ipca_prev = IPCA[433:444]

MSE = 1/length(ipca_prev)*(sum((ipca_prev - f)^2));MSE

MAPE = ((sum(abs((ipca_prev - f)/(ipca_prev))))/length(ipca_prev))*100;MAPE

SMAPE = ((sum(abs((ipca_prev - f))/((ipca_prev + f)/2)))/length(ipca_prev))*100 ;SMAPE

U =sqrt( (sum((ipca_prev - f)^2))/(sum((ipca_prev[2:12] - ipca_prev[1:11])^2)) );U


# Data frame com resultados das métricas
resultado = data.frame(MSE,MAPE,U,SMAPE)
# Gerar uma tabela que será visualizada no documento. Mais detalhes em http://haozhu233.github.io/kableExtra/
knitr::kable(resultado, format = "pandoc", digits = c(2,2,2), align = 'c')


#############################################################################################
#					Histograma das posteriores
#############################################################################################

x11(5.91,5.91)
par(mfrow=c(4,4))
for (i in 1:12)
{
	if (i == 1) {
	hist(gibbsOut.2$dV[cadeia],main="",xlab="V")
	hist(gibbsOut.2$dW[seq(burn,MCMC,thin),i],main="",xlab=paste0("W ",i))
	}else{
	hist(gibbsOut.2$dW[seq(burn,MCMC,thin),i],main="",xlab=paste0("W ",i))
	}
}


x11(6,7)
par(mfrow=c(4,4))
for (i in 1:12)
{
	if (i == 1){
	plot(gibbsOut.2$dV[cadeia],main="",xlab="V",type="l",ylab="")
	plot(gibbsOut.2$dW[cadeia,i],type="l",main="",xlab=paste0("W ",i),ylab="")
	
	}else{
	plot(gibbsOut.2$dW[cadeia,i],type="l",main="",xlab=paste0("W ",i),ylab="")
	hist(gibbsOut.2$dW[seq(burn,MCMC,thin),i],main="",xlab=paste0("W ",i))	
	}
}

mean(gibbsOut.2$dV[cadeia])
plot(gibbsOut.2$dV,type="l")

hist(gibbsOut.2$dV[cadeia],main="",xlab="V")
plot(gibbsOut.2$dV[cadeia],main="",xlab="V",type="l",ylab="")

x11(6,7)
par(mfrow=c(6,2))
for(i in 1:1)
{
	if (i == 1){
	par(mar=c(8,4,8,1))
	hist(gibbsOut.2$dV[cadeia],main="",xlab="V")
	#plot(gibbsOut.2$dV[cadeia],main="",xlab="V",type="l",ylab="")
	#plot(gibbsOut.2$dW[cadeia,i],type="l",main="",xlab=paste0("W ",i),ylab="")
	#hist(gibbsOut.2$dW[seq(burn,MCMC,thin),i],main="",xlab=paste0("W ",i))
	}else{
	plot(gibbsOut.2$dW[cadeia,i],type="l",main="",xlab=paste0("W ",i),ylab="")
	hist(gibbsOut.2$dW[seq(burn,MCMC,thin),i],main="",xlab=paste0("W ",i))	
			}
}

windows(width = 5.91, height = 2,rescale = "fixed")
par(mfrow=c(1,2))
for(i in 1:1)
{

if (i == 1){
	par(mar=c(4.5,4,1.5,4))
	hist(gibbsOut.2$dV[cadeia],main="",xlab=expression(nu^2),ylab="Frequência")
	par(mar=c(4.5,0,1.5,1))
	plot(gibbsOut.2$dV[cadeia],main="",xlab=expression(nu^2),type="l",ylab="")
	par(mar=c(4.5,4,1.5,4))
	hist(gibbsOut.2$dW[seq(burn,MCMC,thin),i],main="",xlab=(expression(omega^2[i])))
	par(mar=c(4.5,0,1.5,1))
	plot(gibbsOut.2$dW[cadeia,i],type="l",main="",xlab=(expression(omega^2[i])),ylab="")
	} else{

	par(mar=c(2,4,1.5,4))
	hist(gibbsOut.2$dW[seq(burn,MCMC,thin),i],main="",xlab=paste0("W ",i))
	par(mar=c(2,0,1.5,1))
	plot(gibbsOut.2$dW[cadeia,i],type="l",main="",xlab=paste0("W ",i),ylab="")
	
	}
}


windows(width = 5.91, height = 2,rescale = "fixed")
par(mfrow=c(1,2))
par(mar=c(4.5,4,1.5,4))
hist(gibbsOut.2$dW[cadeia,12],main="",xlab=(expression(omega[12]^2)),ylab="Frequência")
par(mar=c(4.5,0,1.5,1))
plot(gibbsOut.2$dW[cadeia,12],type="l",main="",xlab=(expression(omega[12]^2)),ylab="Frequência")


#############################################################################################
#					Comparações
#############################################################################################



#############################################################################################
#					Gráfico das distribuições das variâncias
#############################################################################################

x11()
plot(ts.ipca, type="p", ylab="Expenditure - mcmc")
lines(output.2$m[,1],type="l", col = "red")
data<-gibbsOut.2$dW[seq(burn,MCMC,thin),]
teste<-matrix(data,nrow =900 , ncol = 12 )
dim(teste)

w1  <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),1])
w2  <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),2])
w3  <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),3])
w4  <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),4])
w5  <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),5])
w6  <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),6])
w7  <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),7])
w8  <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),8])
w9  <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),9])
w10 <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),10])
w11 <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),11])
w12 <- c(gibbsOut.2$dW[seq(burn,MCMC,thin),12])

data_post <- data.frame(w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12)

q1  <- ggplot(data = data_post , aes(w1)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))
q2  <- ggplot(data = data_post , aes(w2)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))	
q3  <- ggplot(data = data_post , aes(w3)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))
q4  <- ggplot(data = data_post , aes(w4)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))
q5  <- ggplot(data = data_post , aes(w5)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))
q6  <- ggplot(data = data_post , aes(w6)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))
q7  <- ggplot(data = data_post , aes(w7)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))
q8  <- ggplot(data = data_post , aes(w8)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))
q9  <- ggplot(data = data_post , aes(w9)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))
q10 <- ggplot(data = data_post , aes(w10)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))
q11 <- ggplot(data = data_post , aes(w11)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))
q12 <- ggplot(data = data_post , aes(w12)) + geom_histogram() + theme_bw()+
	theme(plot.title=element_text(hjust=0.5))

x11(40,25)	
grid.arrange(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12)


#############################################################################################
#					Gráfico com ajuste e previsão
#############################################################################################

Time = seq(1995,2016,length.out=(2016-1995)*12)
x<-Time
y<-ipca
z<-output.2$m[2:253,1]
ipca2 <- data.frame(x,y,z)

Time_prev <- seq(2016,2017,length.out=(2017 - 2016)*12)
predit<-f
obs <- IPCA[253:264]
data_prev <- data.frame(Time_prev,predit,obs)

grafico <-  ggplot() + geom_line(data = ipca2, aes(x=x,y=z)) +
		geom_point(data = ipca2, aes(x=x,y=y))+
		#geom_ribbon(aes(x = interval_time,ymin = ,ymax = ), fill = "grey70")+
		geom_line(data = data_prev, aes(x = Time_prev,y = Series.1),color = c("blue"))+
		geom_point(data = data_prev, aes(x=Time_prev,y=obs)) +
		theme_bw()+
		theme(plot.title=element_text(hjust=0.5))+
		geom_vline(xintercept = 2016)+
		labs(y=expression("IPCA"), x="Anos") +
		ggtitle(expression(""))

x11(16.18,6.70)
grafico
dev.print(device=pdf, file="DLM_ajuste_prev.pdf")
dev.off()







