rm(list=ls(all=TRUE))

library("forecast")
library("tseries")

##############################################################################
#                         Leitura dos dados do IPCA
##############################################################################

dados<-read.table("IPCA.txt",header=T,sep=",")
dados[,2]
attach(dados)
names(dados)


##############################################################################
#             Análise para o periodo 1995 - 2014 (previsão 12 passos)
##############################################################################

out.sample = 12
ipca<-c(IPCA[181:444])
N = length(ipca)
ts.ipca = ts(ipca[1:(N-out.sample)],frequency=12,start=c(1995))
n = length(ts.ipca)

ipca_fit_arma = arima(ts.ipca, order=c(1,0,12))  
residuo = ipca_fit_arma$residuals
ajustado = ts.ipca - residuo


#previsão
prev<-predict(ipca_fit_arma,out.sample)
ts.ipca.prev = ts(ipca[(N-out.sample+1):N],frequency=12)
real = ts.ipca.prev
F = prev$pred[1:out.sample]

######################### Tabela para métricas Ajuste #######################

MSE = 1/n*(sum((ts.ipca - ajustado)^2));MSE
MAPE = ((sum(abs((ts.ipca - ajustado)/(ts.ipca))))/n)*100;MAPE
SMAPE = ((sum(abs((ts.ipca - ajustado))/((ts.ipca + ajustado)/2)))/n)*100;SMAPE
U =sqrt( (sum((ts.ipca - ajustado)^2))/(sum((ts.ipca[2:n] - ts.ipca[1:(n-1)])^2)) );U
MAD = (sum(abs(ts.ipca - ajustado)))/(n);MAD
resultado = data.frame(MSE,MAPE,SMAPE,MAD,U)
knitr::kable(resultado, format = "pandoc", digits = c(3,3,3), align = 'c',caption="Ajuste")

############### Tabela para métricas previsão para 12 passos #################

MSE = 1/length(real)*(sum((real - F)^2));MSE
MAPE = ((sum(abs((real - F)/(real))))/length(real))*100;MAPE
SMAPE = ((sum(abs((real - F))/((real + F)/2)))/length(real))*100 ;SMAPE
U =sqrt( (sum((real - F)^2))/(sum((real[2:out.sample] - real[1:(out.sample-1)])^2)) );U
MAD = (sum(abs(F - real)))/(length(real));MAD
resultado_prev = data.frame(MSE,MAPE,U,SMAPE)
knitr::kable(resultado_prev, format = "pandoc", digits = c(2,2,2), align = 'c',caption="Previsão")


