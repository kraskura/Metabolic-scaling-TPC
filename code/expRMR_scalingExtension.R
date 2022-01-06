

# 1. Multiple regression model, MTE (Downs et al 2008 Funct Ecol, https://www.jstor.org/stable/20142797)
# ln(MR) = ln(a) + b*ln(BW) + c(1000/T)
# c = -1E/1000k # k = Boltz const. 


k<-(8.62*10^(-5)) # Boltzmann's constant 
# kBoltz<-1.38 * 10^23 
E<-0.63 # activation energy MTE 


# interaction model 
# ln(MR) = ln(a) + b*ln(BW) + c (1000/T) + d (ln(BW) * 1000/T)
# exp(ln(MR)) = exp(ln(a) + b*ln(BW) + c (1000/T) + d (ln(BW) * 1000/T)

# ----------------------------------------------------------------------------------------------------------------

# RELATEABLE WORK, builds up
# Metabolic scaling metadata analyses multiregression models: 


# RMR_MultReg_ER_model4 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
# RMR_MultReg_ER_model4.o <- lmer(lnRMR ~ offset(1 * lnBWg) + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
# RMR_MultReg_ER_model4int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
# 
# dataPred.model4<-expand.grid(lnBWg = seq(log(1), log(100)), tempTest = seq(0,35, 5))
# dataPred.model4$tempTestK<-celsius.to.kelvin(dataPred.model4$tempTest, round = 2)
# dataPred.model4$tempTestK1<-1/dataPred.model4$tempTestK
# dataPred.model4$tempTestK1000<-1000/dataPred.model4$tempTestK
# 
# dataPred.model4$pred <- predict(RMR_MultReg_ER_model4, newdata = dataPred.model4, re.form = NA)
# dataPred.model4$pred.o <- predict(RMR_MultReg_ER_model4.o, newdata = dataPred.model4, re.form = NA)
# dataPred.model4$pred.i <- predict(RMR_MultReg_ER_model4int, newdata = dataPred.model4, re.form = NA)
# 
# # backtransform:
# dataPred.model4$MR<-exp(dataPred.model4$pred)
# dataPred.model4$MR.o<-exp(dataPred.model4$pred.o)
# dataPred.model4$MR.i<-exp(dataPred.model4$pred.i)
# 
# dataPred.model4$BWg<-exp(dataPred.model4$lnBWg)
# 
# ggplot(dataPred.model4, aes(x = tempTest, y = MR.o, size = BWg, group = BWg))+
#   geom_line(mapping = aes(x = tempTest, y = MR, size = BWg, group = BWg), color = "blue", size = 0.5, lty=2)+
#   geom_line(mapping = aes(x = tempTest, y = MR.i, size = BWg, group = BWg), color = "red", size = 1, lty = 3)+
#   geom_line(size=1)+
#   theme_light()

# ----------------------------------------------------------------------------------------------------------------


# Applying this model structure here. This can be easily used by reserachers and would be accessibel format. 

# 1. Simple model
# ln(MR) = ln(a) + b*ln(BW) + c (1000/T) 

# 2. Interaction model (all same idea, rewritten in different formats)
# ln(MR) = ln(a) + b*ln(BW) + c (1000/T) + d (ln(BW) * 1000/T)
# ln(MR)  = ln(a) + b*ln(BW) + c*(1000/T) + d*(ln(BW) * 1000/T)
# ln(MR)  = ln(a) + (b + d*(1000/T))*ln(BW) + c*(1000/T)
# ln(MR)  = ln(a) + b*ln(BW) + (c + d*ln(BW))*(1000/T)

# Some key parameters 
a <- 17
b <- 0.89
E <- 0.63
  
# create 
dd<-expand.grid(BW=seq(100,1000, 100), K = (celsius.to.kelvin(seq(5,25))), d=c(-0.1, -0.01, -0.05), E = c(0.5, 0.63, 0.2))
dd$lnMR1<-0
dd$lnMRi<-0

# model with scaling and temperature effect using Arrh activ energy (core MTE model)
simple.exp.mod1<-function(a, b, E, k_Boltz = (8.62*10^(-5)), K, BW){
  c<-(-1*E)/(1000*k_Boltz)
  lnMR<- a + (b*log(BW)) + (c*(1000/K))
  return(lnMR)
}

# model with scaling and temperature effect, but interaction betweem temp and mass
# the model would suggest that the rate of exponential increase of MR with temp (e.g. accelerated increase) will be greater with increasing body mass. 
i.exp.mod1<-function(a, b, E, k_Boltz = (8.62*10^(-5)), K, BW, d){
  c<-(-1*E)/(1000*k_Boltz)
  # lnMR<- a + (b*log(BW)) + (c*(1000/K)) + (d*((1000/K) * (log(BW))))
  lnMR  = a + (b*log(BW)) + ((c + (d*log(BW))) * (1000/K)) # same function, re-written. 
  return(lnMR)
}

for(i in 1:nrow(dd)){
  dd$lnMR1[i]<-simple.exp.mod1(a=15, b=0.89, E=0.6, K = dd$K[i], BW = dd$BW[i])
  dd$lnMRi[i]<-i.exp.mod1(a=15, b=0.75, E=dd$E[i], K = dd$K[i], BW = dd$BW[i], d =dd$d[i])
}


# simple model plotted , no interaction
ggplot(dd, aes(x = K, y = exp(lnMR1)/BW, size = BW, group = BW))+
  geom_line(size=1)+
  theme_light()
  
# interaction model, plotted mass specific 
ggplot(dd, aes(x = K, y = exp(lnMRi)/BW, group = interaction(E,d), color = d))+
  geom_line(size=1)+
  theme_light()+
  facet_wrap(.~BW, scales = "free")

# top numbers in each panel are E (acrtiv energy.) and bottom is d (the interaction term in model)
ggplot(dd, aes(x = K, y = exp(lnMRi)/BW, group = BW, color = BW))+
  geom_line(size=1)+
  theme_light()+
  facet_wrap(E~d, scales = "free")








