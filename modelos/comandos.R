rm(list=ls(all=TRUE))

library("forecast")
library("tseries")
#library("ggplot2")
#library("gridExtra")
#library("grDevices")
#library("grid")
#library(tidyverse)
#library("MASS")
#library("nortest")
#library("car")

####ipca minusculo = serie de 1995 ate 2015
####IPCA maiusculo = serie completa.
####IPCA[181:444] = 1995 ate 2016
####IPCA[181:420] 
##Dados IPCA

##############################################################################
#                         Leitura dos dados do IPCA
##############################################################################

dados<-read.table("IPCA.txt",header=T,sep=",")
dados[,2]
attach(dados)
names(dados)


##############################################################################
#             Dados corrigidos para o periodo correto 1995 - 2015
##############################################################################

ipca<-c(IPCA[181:432])
ts.ipca=ts(ipca,frequency=12,start=c(1995)) #start=c(1994,7)
 
##############################################################################
#                     Gráfico de IPCA 1995 - 2016
##############################################################################


#x11(16.06,9.79)
x11(3.94,3.15)
par(mar=c(4,4,1,1))
plot(IPCA[181:444],main="",type="l",xlab="",,xlim=c(0,276),ylab="",axes=FALSE,lwd=2)#IPCA 1994 - 2016
title(ylab="IPCA - Variação Mensal",cex.lab=1)
title(xlab="Período",cex.lab=1)
#title(main="IPCA no periodo de 1995 - 2016",cex.main=2)
axis(1,seq(1,276,12),(1995:2017),cex.axis=1)
axis(2,cex.axis=1)

			########### Gráfico com GGPLOT ###########

Time = seq(1995,2016,length.out=(2016-1995)*12)
x<-Time
y<-ipca
ipca2 <- data.frame(x,y)
ipca2
q3<-  ggplot()+ geom_line(data = ipca2, aes(x=x,y=y)) +
	theme_bw()+
	theme(plot.title=element_text(hjust=0.5))+
	labs(y=expression("IPCA"), x="Anos") +
	ggtitle(expression(""))
x11(16.06,5.79)
q3
dev.print(device=pdf, file="ipca.pdf")
dev.off()

##############################################################################
#              Gráficos de Autocorrelação e Autocorrelação Parcial
##############################################################################

x11(16.188976,6.69291)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
#x11(5.75,4.5)
acf(ipca,lag.max=30,ylab="FAC",main="",xlab="lag",axes=FALSE)
#title(main="Função de autocorrelação",cex.main=1.5)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)
pacf(ipca,lag.max=30,main="",cex.main=2,ylab="FACP",xlab="lag",axes=FALSE)
#title(main="Função de autocorrelação parcial",cex.main=1.5)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)

			########### Gráfico com GGPLOT ###########

Time<-seq(0,29,1)
bacf<-acf(ipca,lag.max=30)
bpacf<-pacf(ipca,lag.max=30)
bacfdf <- with(bacf, data.frame(lag, acf))
bpacfdf <- with(bpacf, data.frame(lag, acf))

q1 <- #ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
      ggAcf(ipca, lag.max=30) +
	#geom_hline(aes(yintercept = 0)) +
      #geom_segment(mapping = aes(xend = lag, yend = 0))+
	theme_bw()+
	theme(plot.title=element_text(hjust=0.5))+
	labs(y=expression("IPCA"), x="Lags") +
	ggtitle(expression("Gráfico de função de autocorrelação"))

q2 <- #ggplot(data = bpacfdf , mapping = aes(x = lag, y = acf)) +
      ggPacf(ipca, lag.max=30)+
	#geom_hline(aes(yintercept = 0)) +
      #geom_segment(mapping = aes(xend = lag, yend = 0))+
	theme_bw()+
	theme(plot.title=element_text(hjust=0.5))+
	labs(y=expression("IPCA"), x="Lags") +
	ggtitle(expression("Gráfico de função de autocorrelação parcial"))

grid.arrange(q1,q2)
dev.print(device=pdf, file="autocorrelation.pdf")
dev.off()

## Analisando os correlogramas. Sugestão de Modelo ARMA: 
## q=14 ou 5 ou 4(medias moveis); p= 1(autoregressiva), 
## Pela decomposição da série: tem sazonalidade anual

##############################################################################
#                         Testes de Estacionariedade
##############################################################################

############## (se precisasse estacionarizar, diff((ipca)) ###################

adf.test(ts.ipca) # Dickey-Fuller 
pp.test(ts.ipca) # Phillips-Perron

##############################################################################
#                         Decomposição da série
##############################################################################


x11(5.11,5.11)
#x11(4.07,4.29)
par(mar=c(1,2,1,2))
plot(stl((ts.ipca), "per"),labels=rbind("Observada","Sazonalidade","Tendência","Ruído"))
plot(stl((ts.ipca), "per"),main="teste",labels=cbind("teste"))


decomposta=stl(ts.ipca, "per")
#decomposta

decompose.ipca=decompose(ts.ipca)
plot(decompose.ipca)
names(decompose.ipca)
decompose.ipca$x
decompose.ipca$seasonal
decompose.ipca$trend
decompose.ipca$random

##############################################################################
#                         Ajustando o Modelo ARMA
##############################################################################


ipca_fit_arma=arima(ts.ipca, order=c(1,0,12)) # q=4, 5 ou 14, period=12  

ipca_fit_arma$coef
ipca_fit_arma$sigma2  #510810
confint(ipca_fit_arma)
residuo=ipca_fit_arma$residuals


acf(residuo,lag.max=30,main="",cex.main=2,ylab="FACP",xlab="lag",axes=FALSE)
title(main="Função de Autocorrelação dos Residuos - ARMA(1,12)",cex.main=1)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)

ipca_fit_arma4=arima(ts.ipca, order=c(1,0,4)) # q=4, 5 ou 14, period=12
residuo4 = ipca_fit_arma4$residuals
ipca_fit_arma5=arima(ts.ipca, order=c(1,0,5)) # q=4, 5 ou 14, period=12
residuo5 = ipca_fit_arma5$residuals

x11(16.188976,6.69291)
par(mfrow=c(1,2))
x11(5.75,3)
acf(residuo,lag.max=30,ylab="FAC",main="",xlab="lag",axes=FALSE)
title(main="Função de Autocorrelação dos Residuos - ARMA(1,12)",cex.main=1)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)
pacf(residuo,lag.max=30,main="",cex.main=2,ylab="FACP",xlab="lag",axes=FALSE)
title(main="Função de Autocorrelação Parcial dos Residuos - ARMA(1,12)",cex.main=1)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)


ipca_fit_arma1=arima(ts.ipca, order=c(1,0,1)) # q=4, 5 ou 14, period=12
residuo1 = ipca_fit_arma1$residuals

acf(residuo1,lag.max=30,main="",cex.main=2,ylab="FACP",xlab="lag",axes=FALSE)
title(main="Função de Autocorrelação dos Residuos - ARMA(1,1)",cex.main=1)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)

##############################################################################
#                         		Previsões
##############################################################################

#previsão para 12 meses
prev<-predict(ipca_fit_arma,12)


##############################################################################
#                         Critérios de Comparação de Modelos
##############################################################################

ipca_fit_arma=arima(ts.ipca, order=c(1,0,12)) # q=4, 5 ou 14, period=12  
residuo=ipca_fit_arma$residuals

jarque.bera.test(residuo)


aic=ipca_fit_arma$aic
#mape=((1/length(ts.ipca))*sum(abs(residuo)/ts.ipca)*100)
smape=100*sum((abs(residuo))/((2*ts.ipca-residuo)/2))*(1/length(ts.ipca))#coloquei a subtração
utheil=sqrt((sum(residuo)^2))/(sqrt(sum(ts.ipca)+sqrt(sum(ts.ipca+residuo)^2))) #verificar esta equação
print(rbind(aic,smape,utheil))

real = IPCA[433:444]
residuo1 = prev$pred-IPCA[253:264]
ipca_prev = real
F = prev$pred
ajustado = ipca-residuo

mad = (1/length(prev$pred))*sum(abs(prev$pred-IPCA[253:264]))
mse = (1/length(prev$pred))*sum((prev$pred-IPCA[253:264])^2)
mape = 100*((1/length(prev$pred))*sum(abs(residuo1/real)))

MSE = 1/length(ts.ipca)*(sum((ts.ipca - ajustado)^2));MSE

MAPE = ((sum(abs((ts.ipca - ajustado)/(ts.ipca))))/length(ts.ipca))*100;MAPE

SMAPE = ((sum(abs((ts.ipca - ajustado))/((ts.ipca + ajustado)/2)))/length(ts.ipca))*100

U =sqrt( (sum((ts.ipca - ajustado)^2))/(sum((ts.ipca[2:252] - ts.ipca[1:251])^2)) );U

MAD = (sum(abs(ts.ipca - ajustado)))/(length(ts.ipca));MAD

# Data frame com resultados das métricas
resultado = data.frame(MSE,MAPE,U,SMAPE)
# Gerar uma tabela que será visualizada no documento. Mais detalhes em http://haozhu233.github.io/kableExtra/
knitr::kable(resultado, format = "pandoc", digits = c(2,2,2), align = 'c')



##############################################################################
MSE = 1/length(real)*(sum((real - F)^2));MSE

MAPE = ((sum(abs((ipca_prev - F)/(ipca_prev))))/length(ipca_prev))*100;MAPE

SMAPE = ((sum(abs((ipca_prev - F))/((ipca_prev + F)/2)))/length(ipca_prev))*100 ;SMAPE


U =sqrt( (sum((ipca_prev - F)^2))/(sum((ipca_prev[2:12] - ipca_prev[1:11])^2)) );U

MAD = (sum(abs(F - real)))/(length(real));MAD


# Data frame com resultados das métricas
resultado = data.frame(MSE,MAPE,U,SMAPE)
# Gerar uma tabela que será visualizada no documento. Mais detalhes em http://haozhu233.github.io/kableExtra/
knitr::kable(resultado, format = "pandoc", digits = c(2,2,2), align = 'c')




bic=AIC(ipca_fit_arma,k = log(length(ts.ipca)))
bic
#12 > bic
[1] 204.3133
> bic 
[1] 204.3133
> aicc 
[1] 153.4058
smape > smape
[1] 34.74624

#4
> bic
[1] 186.0612
> bic 
[1] 186.0612
> aicc 
[1] 161.8142
44.82196199
#5> bic
[1] 187.8595
> bic 
[1] 187.8595
> aicc 
[1] 160.2166

npar <- length(ipca_fit_arma$coef) + 1
nstar <- length(ipca_fit_arma$residuals) - ipca_fit_arma$arma[6] - ipca_fit_arma$arma[7] * ipca_fit_arma$arma[5]

bic <- ipca_fit_arma$aic + npar * (log(nstar) - 2)
aicc <- ipca_fit_arma$aic + 2 * npar * (nstar/(nstar - npar - 1) - 1)



#opção para o utheil
n=length(ts.ipca)
utheil_novo=sqrt((sum((residuo[2:n])^2))/(sum((ts.ipca[2:n]-ts.ipca[1:n-1])^2)))

utheil_novo=sqrt((sum((residuo[2:n])^2))/(sum((-residuo[2:n]+ts.ipca[2:n]-ts.ipca[1:n-1])^2)))




#gráfico com os valores ajustados.

ajustado<-(ipca-residuo)
ajustado
x<-seq(1,length(ipca),1)
y<-1:264
k<-253:264

x11(5.75,3.54)
plot(IPCA[181:444],main="",type="l",xlab="",ylab="",axes=FALSE,lwd=2)#IPCA 1994 - 2016
title(ylab="IPCA Variação Mensal",cex.lab=1.5)
title(xlab="Período",cex.lab=1.5)
axis(1,seq(1,length(IPCA[181:444]),12),(1995:2016),lwd=1.6,cex.axis=1.4)
axis(2,lwd=1.6,cex.axis=1.4)

##############################################################################
#		GRÁFICO DE AJUSTE COM OBSERVADO E PREVISTO
##############################################################################

#x11(16.18,6.70)
windows(width = 5.9, height = 4, rescale = "R")
par(mar=c(4,4,2,0))#bottom,left,top,rigth
plot(ipca,pch=3,axes=FALSE,ylab="",xlab="",xlim=c(0,276),cex.lab=1,cex=1,cex.axis=1)
axis(1,at=seq(1,276,12),(1995:2017),cex.axis=1)
axis(2,cex.axis=1)
title(ylab="IPCA",cex.lab=1)
title(xlab="Período",cex.lab=1)
#title(main="Ajuste e previsão do IPCA pelo ARIMA(1,0,12)",cex.lab=1)
lines(x,ajustado,type="l",lty=1,lwd=2)
lines(k,IPCA[433:444],type="p",cex.lab=1.4,pch=19)
abline(v=252,lty=2)
lines(k,prev$pred,lty=1,lwd=2,col="blue")
lines(k,prev$pred+1.96*prev$se,lty=2)
lines(k,prev$pred-1.96*prev$se,lty=2)

		################ GRÁFICO COM GGPLOT ################ 



Time_prev <- seq(2016,2017,length.out=(2017 - 2016)*12)
Time_prev
predit<-prev$pred
obs <- IPCA[253:264]
data_prev <- data.frame(Time_prev,predit,obs)

lwr <- prev$pred-1.96*prev$se
upr <- prev$pred+1.96*prev$se
interval_time<- Time_prev
interval<-data.frame(lwr,upr,interval_time)

q3<-  ggplot()+ geom_line(data = ipca2, aes(x=x,y=z)) +
	geom_point(data = ipca2, aes(x=x,y=y)) +
	geom_ribbon(aes(x = interval_time,ymin= lwr,ymax=upr),fill = "grey70" ) +
	geom_line(data = data_prev, aes(x=Time_prev,y=predit),color = c("blue")) + 
	geom_point(data = data_prev, aes(x=Time_prev,y=obs)) +
	theme_bw()+
	theme(plot.title=element_text(hjust=0.5))+
	geom_vline(xintercept = 2016) +
	labs(y=expression("IPCA"), x="Anos") +
	ggtitle(expression(""))
x11(16.18,6.70)
q3
dev.print(device=pdf, file="Arma_ajuste_prev.pdf")
dev.off()
#
#				Parar aqui....
##############################################################################
plot(IPCA[253:264],type="l",cex.lab=1.4,pch=19)
lines(1:12,prev$pred,lty=3,lwd=2)

#Felipe, faça esses gráficos todos e salve organizando com nomes. Comece escrever um documento que contenha os dados e seus gráficos, explicando os dados e a finalidade do trabalho.
#Aplique esses testes todos e faça tabelas com os resultados.


########################## ARMA ATE 2014 #################

###GRÁFICOS DE AUTOCORRELAÇÃO###

par(mfrow=c(2,1))
acf(IPCA[181:420],lag.max=30,ylab="FAC",main="",xlab="lag",axes=FALSE)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)
title(main="Correlograma de 1995 até 2014")
pacf(IPCA[181:420],lag.max=30,main="",ylab="FACP",xlab="lag",axes=FALSE)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)


#Criando objeto ts e decompondo a série 
ts.IPCA=ts(IPCA[181:420],frequency=12,start=c(1995)) #start=c(1994,7)

#Testes de Estacionariedade (se precisasse estacionarizar, diff((ipca))
adf.test(ts.IPCA) # Dickey-Fuller 

#Ajustando o Modelo ARMA
IPCA_fit_arma=arima(ts.IPCA, order=c(1,0,12)) # q=4, 5 ou 14, period=12  

IPCA_fit_arma$coef
IPCA_fit_arma$sigma2  #510810
confint(IPCA_fit_arma)
Residuo=IPCA_fit_arma$residuals

#Critérios de Comparação de Modelos
aic=IPCA_fit_arma$aic
#mape=((1/length(ts.IPCA))*sum(abs(Residuo)/ts.IPCA)*100)
smape=100*sum((abs(Residuo))/((2*ts.IPCA-Residuo)/2))*(1/length(ts.IPCA))#coloquei a subtração
utheil=sqrt((sum(Residuo)^2))/(sqrt(sum(ts.IPCA)+sqrt(sum(ts.IPCA+Residuo)^2))) #verificar esta equação
print(rbind(aic,smape,utheil))

PREV<-predict(IPCA_fit_arma,24)
PREV

###GRÁFICO DO AJUSTE E PREVISÃO###

ajustado<-(IPCA[181:420]-Residuo)
ajustado
x<-seq(1,length(IPCA[181:420]),1)
y<-1:264
k<-241:264
plot(IPCA[181:420],pch=3,axes=FALSE,ylab="",xlab="",xlim=c(0,276))
axis(1,at=seq(1,276,12),(1995:2017),lwd=2,cex.lab=1.4,cex.axis=1.4)
axis(2,lwd=2,cex.axis=1.4,cex.axis=1.4)
title(ylab="IPCA Variação Mensal",cex.lab=1.4)
title(xlab="Período",cex.lab=1.4)
lines(x,ajustado,type="l",cex.lab=1.4)
lines(k,IPCA[241:264],type="p",cex.lab=1.4,pch=19)
abline(v=240,lty=2)
lines(k,PREV$pred,lty=1,lwd=2,col="blue")
lines(k,PREV$pred+1.96*PREV$se,lty=2)
lines(k,PREV$pred-1.96*PREV$se,lty=2)
title(main="ARMA(1,0,12) de 1995 ate 2014")

########################## ARMA ATE 2012 #################

###GRÁFICOS DE AUTOCORRELAÇÃO###

ipca_dif<-diff(IPCA[181:396])

par(mfrow=c(1,2))
acf(ipca_dif,lag.max=30,ylab="FAC",main="",xlab="lag",axes=FALSE)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)
pacf(ipca_dif,lag.max=30,main="",ylab="FACP",xlab="lag",axes=FALSE)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)


#Criando objeto ts e decompondo a série 
ts.IPCA=ts(ipca_dif,frequency=12,start=c(1995)) #start=c(1994,7)

#Testes de Estacionariedade (se precisasse estacionarizar, diff((ipca))
adf.test(ipca_dif) # Dickey-Fuller 

#Ajustando o Modelo ARMA
IPCA_fit_arma=arima(ts.IPCA, order=c(2,0,2)) # q=4, 5 ou 14, period=12  

IPCA_fit_arma$coef
IPCA_fit_arma$sigma2  #510810
confint(IPCA_fit_arma)
Residuo=IPCA_fit_arma$residuals

#Critérios de Comparação de Modelos
aic=IPCA_fit_arma$aic
#mape=((1/length(ts.IPCA))*sum(abs(Residuo)/ts.IPCA)*100)
smape=100*sum((abs(Residuo))/((2*ts.IPCA-Residuo)/2))*(1/length(ts.IPCA))#coloquei a subtração
utheil=sqrt((sum(Residuo)^2))/(sqrt(sum(ts.IPCA)+sqrt(sum(ts.IPCA+Residuo)^2))) #verificar esta equação
print(rbind(aic,smape,utheil))

PREV<-predict(IPCA_fit_arma,48)
PREV

###GRÁFICO DO AJUSTE E PREVISÃO###

ajustado<-(ipca_dif-Residuo)
ajustado
x<-seq(1,length(ipca_dif),1)
y<-1:264
k<-217:264
plot(IPCA[181:396],pch=3,axes=FALSE,ylab="",xlab="",xlim=c(0,276))
axis(1,at=seq(1,276,12),(1995:2017),lwd=2,cex.lab=1.4,cex.axis=1.4)
axis(2,lwd=2,cex.axis=1.4,cex.axis=1.4)
title(ylab="IPCA Variação Mensal",cex.lab=1.4)
title(xlab="Período",cex.lab=1.4)
lines(x,ajustado,type="l",cex.lab=1.4)
lines(k,IPCA[217:264],type="l",cex.lab=1.4,pch=19)
abline(v=216,lty=2)
lines(k,PREV$pred,lty=3,lwd=2)
lines(k,PREV$pred+1.96*PREV$se,lty=2)
lines(k,PREV$pred-1.96*PREV$se,lty=2)







par(mar=c(2,4,2,2),cex.lab=1.4,cex=1.4,cex.main=1.4,cex.axis=1.4)
layout( matrix(c(1,1,2,2,3,3,4,4),4,2,byrow=TRUE),c(1,1,1,1),c(1,1,1,1),respect=TRUE)

#mar c(bottom, left, top, right) 
par(mar=c(3,0,1,0),cex.lab=1.3,cex=1.2,cex.main=1.2,cex.axis=1)
plot(ipca,main="",type="l",xlab="",ylab="",axes=FALSE,lwd=1)#IPCA 1994 - 2016
title(ylab="IPCA Variação Mensal",cex.lab=1.5)
#title(xlab="Período",cex.lab=1.5)
#axis(1,seq(1,length(ipca),12),(1994:2016),lwd=1.6,cex.axis=1.4)
axis(2,lwd=1.6,cex.axis=1)

par(mar=c(3,0,1,0),cex.lab=1.3,cex=1.2,cex.main=1.2,cex.axis=1)
plot(decompose.ipca$seasonal,main="",type="l",xlab="",ylab="",axes=FALSE,lwd=1)#IPCA 1994 - 2016
title(ylab="IPCA Variação Mensal",cex.lab=1.5)
#title(xlab="Período",cex.lab=1.5)
#axis(1,seq(1,length(ipca),12),(1994:2016),lwd=1.6,cex.axis=1.4)
axis(2,lwd=1.6,cex.axis=1)

par(mar=c(3,0,1,0),cex.lab=1.3,cex=1.2,cex.main=1.2,cex.axis=1)
plot(decompose.ipca$trend,main="",type="l",xlab="",ylab="",axes=FALSE,lwd=1)#IPCA 1994 - 2016
title(ylab="IPCA Variação Mensal",cex.lab=1.5)
#title(xlab="Período",cex.lab=1.5)
#axis(1,seq(1,length(ipca),12),(1994:2016),lwd=1.6,cex.axis=1.4)
axis(2,lwd=1.6,cex.axis=1)

par(mar=c(3,0,1,0),cex.lab=1.3,cex=1.2,cex.main=1.2,cex.axis=1)
plot(decompose.ipca$random,main="",type="l",xlab="",ylab="",axes=FALSE,lwd=1)#IPCA 1994 - 2016
title(ylab="IPCA Variação Mensal",cex.lab=1.5)
title(xlab="Período",cex.lab=1.5)
axis(1,seq(1,length(ipca),12),(1994:2016),lwd=1.6,cex.axis=1.4)
axis(2,lwd=1.6,cex.axis=1)


decompose.ipca$x
decompose.ipca$seasonal
decompose.ipca$trend
decompose.ipca$random


############## IPCA DE 2000 ATÉ 2012 ###########


###GRÁFICOS DE AUTOCORRELAÇÃO###

ipca_2000<-IPCA[241:396]

par(mfrow=c(2,1))
acf(ipca_2000,lag.max=30,ylab="FAC",main="",xlab="lag",axes=FALSE)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)
title(main="correlogramas para IPCA iniciando em 2000 até 2012")
pacf(ipca_2000,lag.max=30,main="",ylab="FACP",xlab="lag",axes=FALSE)
axis(1,lwd=1,cex.axis=1)
axis(2,lwd=1,cex.axis=1)

#Criando objeto ts e decompondo a série 
ts.IPCA=ts(ipca_2000,frequency=12,start=c(2000)) #start=c(1994,7)

#Testes de Estacionariedade (se precisasse estacionarizar, diff((ipca))
adf.test(ipca_2000) # Dickey-Fuller 
pp.test(ipca_2000)

#Ajustando o Modelo ARMA
IPCA_fit_arma=arima(ts.IPCA, order=c(1,0,12)) # q=4, 5 ou 14, period=12  

IPCA_fit_arma$coef
IPCA_fit_arma$sigma2  #510810
confint(IPCA_fit_arma)
Residuo=IPCA_fit_arma$residuals

#Critérios de Comparação de Modelos
aic=IPCA_fit_arma$aic
#mape=((1/length(ts.IPCA))*sum(abs(Residuo)/ts.IPCA)*100)
smape=100*sum((abs(Residuo))/((2*ts.IPCA-Residuo)/2))*(1/length(ts.IPCA))#coloquei a subtração
utheil=sqrt((sum(Residuo)^2))/(sqrt(sum(ts.IPCA)+sqrt(sum(ts.IPCA+Residuo)^2))) #verificar esta equação
print(rbind(aic,smape,utheil))

n=length(ts.IPCA)
utheil_novo=sqrt((sum((Residuo[2:n])^2))/(sum((ts.IPCA[2:n]-ts.IPCA[1:n-1])^2)))
utheil_novo

PREV<-predict(IPCA_fit_arma,48)
PREV

###GRÁFICO DO AJUSTE E PREVISÃO###

ajustado<-(ipca_2000-Residuo)
ajustado
x<-seq(1,length(ipca_2000),1)
y<-1:264
k<-157:204
plot(ipca_2000,pch=3,axes=FALSE,ylab="",xlab="",xlim=c(0,216))
axis(1,at=seq(1,216,12),(2000:2017),lwd=2,cex.lab=1.4,cex.axis=1.4)
axis(2,lwd=2,cex.axis=1.4,cex.axis=1.4)
title(ylab="IPCA Variação Mensal",cex.lab=1.4)
title(xlab="Período",cex.lab=1.4)
lines(x,ajustado,type="l",cex.lab=1.4)
lines(k,IPCA[217:264],type="p",cex.lab=1.4,pch=19)
abline(v=156,lty=2)
lines(k,PREV$pred,lty=1,lwd=2,col="blue")
lines(k,PREV$pred+1.96*PREV$se,lty=2)
lines(k,PREV$pred-1.96*PREV$se,lty=2)
title(main="ARMA(1,0,4) 2000 ate 2012")

##############################################################################
# Get the time values for the time series
decomposta <- decompose(ts.ipca)
Time = seq(1995,2016,length.out=(2016-1995)*12)

# Convert td to data frame
dat = cbind(Time, with(decomposta, data.frame(Observedo=x, Tendência=trend, Sazonalidade=seasonal, Ruído=random)))

plot<- ggplot(gather(dat, component, value, -Time), aes(Time, value)) +
  facet_grid(component ~ ., scales="free_y") +
  geom_line() +
  theme_bw() +
  labs(y=expression("IPCA"), x="Anos") +
  ggtitle(expression("Decomposição da serie temporal IPCA")) +
  theme(plot.title=element_text(hjust=0.5))

dev.print(device=pdf, file="price.pdf")
dev.off()

##############################################################################


#############################################################################

	