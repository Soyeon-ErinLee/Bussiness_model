data =read.csv("iotdata.csv")
data = data[,-c(1,7,9)]
head(data)
data$floor=as.numeric(as.factor(data$floor))
data$roomid=factor(data$roomid)
attach(data)
summary(data)
round(cor(sapply(data, as.integer)), 3)
corrplot::corrplot(round(cor(sapply(data, as.integer)), 3), tl.col = "black", title = "Pearson Correlations", mar = c(4, 0, 4, 0))
hist(humidity)
par(mfrow=c(1,1))
pie(c(length(which(door==1)),length(which(door==0))),labels=c("1","0"), main = "Pie Chart of Door")
hist(floor, breaks = 3)
pie(c(length(which(floor==1)),length(which(floor==2)),length(which(floor==3)),length(which(floor==4)),length(which(floor==5))),labels=c("floorB","floor2","floor3","floor4","floor5"), main = "Pie Chart of Floor")
pie(c(length(which(SSS==1)),length(which(SSS==2)),length(which(SSS==3)),length(which(SSS==4)),length(which(SSS==5))),labels=c("1","2","3","4","5"), main = "Pie Chart of SSS")

# SSS missing
cor(SSS, temperature, method='spearman') # 0.3763913
shapiro.test(SSS); shapiro.test(temperature) # 둘다 normal 아님
m_sss = 0.376*temperature
m_sss = (m_sss - min(m_sss))/(max(m_sss)-min(m_sss))
set.seed(100)
na.ind = sample(which(m_sss>quantile(m_sss, 0.5)), nrow(data)*0.15)
wilcox.test(temperature, temperature[-na.ind]) #  p-value = 5.449e-05
mis = subset(data, select = c("SSS","temperature"))
mis$SSS[na.ind] = NA
mis$Missing = ifelse(is.na(mis$SSS), "Missing", "Observed")
p1 = ggplot(data = mis) + geom_histogram(aes(x=temperature,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 20) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of temp (15% missing)") # 3개 그래프 뒤에서 다같이 그릴 것

set.seed(100)
na.ind = sample(which(m_sss>quantile(m_sss, 0.5)), nrow(data)*0.3)
wilcox.test(temperature, temperature[-na.ind]) #  p-value < 2.2e-16
mis$SSS[na.ind] = NA
mis$Missing = ifelse(is.na(mis$SSS), "Missing", "Observed")
p2 = ggplot(data = mis) + geom_histogram(aes(x=temperature,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 20) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of temp (30% missing)") # missing 전후 density 비교

set.seed(100)
na.ind = sample(which(m_sss>quantile(m_sss, 0.5)), nrow(data)*0.45)
par(mfrow=c(2,1))
wilcox.test(temperature, temperature[-na.ind]) #  p-value < 2.2e-16
mis$SSS[na.ind] = NA
mis$Missing = ifelse(is.na(mis$SSS), "Missing", "Observed")
p3 = ggplot(data = mis) + geom_histogram(aes(x=temperature,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 20) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of temp (45% missing)") # missing 전후 density 비교
require(gridExtra)
grid.arrange(p1, p2, p3, ncol=3)


# fire missing
p1 = hist(humidity)
p2 = hist(humidity[-na.ind])
plot(p1,col=rgb(0,0,1,1/4))
plot(p2, col=rgb(1,0,0,1/4),add=T)
wilcox.test(humidity, humidity[-na.ind]) # p-value  0.0003115
wilcox.test(temperature, temperature[-na.ind]) # p-value = 0.005157

cor(fire, humidity) # -0.374
cor(fire, temperature) # 0.453
m_sss = -0.374*humidity[-na.ind]+0.453*temperature[-na.ind]   
m_sss = (m_sss - min(m_sss))/(max(m_sss)-min(m_sss))
set.seed(100)
na.ind2 = sample(which(m_sss>quantile(m_sss, 0.5)), nrow(data[-na.ind,])*0.15)
mis = subset(data[-na.ind], select = c("fire","humidity", "temperature"))
mis$fire[-na.ind2] = NA
mis$Missing = ifelse(is.na(mis$fire), "Missing", "Observed")
p1_1 = ggplot(data = mis) + geom_histogram(aes(x=humidity,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 30) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of humidity (15% missing with 5% humidity missing)") # missing 전후 density 비교
plot(p1_1)

set.seed(100)
na.ind2 = sample(which(m_sss>quantile(m_sss, 0.5)), nrow(data[-na.ind,])*0.30)
mis = subset(data[-na.ind], select = c("fire","humidity", "temperature"))
mis$fire[-na.ind2] = NA
mis$Missing = ifelse(is.na(mis$fire), "Missing", "Observed")
p1_2 = ggplot(data = mis) + geom_histogram(aes(x=humidity,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 30) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of humidity (30% missing with 5% humidity missing)") # missing 전후 density 비교


set.seed(100)
na.ind2 = sample(which(m_sss>quantile(m_sss, 0.5)), nrow(data[-na.ind,])*0.45)
mis = subset(data[-na.ind], select = c("fire","humidity", "temperature"))
mis$fire[-na.ind2] = NA
mis$Missing = ifelse(is.na(mis$fire), "Missing", "Observed")
p1_3 = ggplot(data = mis) + geom_histogram(aes(x=humidity,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 30) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of humidity (45% missing with 5% humidity missing)") # missing 전후 density 비교


set.seed(100)
na.ind = sample(which(m_sss>quantile(m_sss, 0.5)), nrow(data)*0.45)
wilcox.test(humidity, humidity[-na.ind]) # p-value < 2.2e-16
wilcox.test(temperature, temperature[-na.ind]) # p-value < 2.2e-16
mis$fire[na.ind] = NA
mis$Missing = ifelse(is.na(mis$fire), "Missing", "Observed")
p1_3=ggplot(data = mis) + geom_histogram(aes(x=humidity,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 20) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of humidity (45% missing)") # missing 전후 density 비교
p2_3=ggplot(data = mis) + geom_histogram(aes(x=temperature,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 20) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of temp (45% missing)") # missing 전후 density 비교

require(gridExtra)
grid.arrange(p1_1, p1_2, p1_3, ncol=3) # humidity

grid.arrange(p2_1, p2_2, p2_3, ncol=3) # temp


# co2 missing
shapiro.test(co2)
cor(co2, dust) # -0.3405261
m_sss = 0.341*dust
m_sss = (m_sss - min(m_sss))/(max(m_sss)-min(m_sss))
set.seed(100)
na.ind = sample(which(m_sss>quantile(m_sss, 0.5)), nrow(data)*0.15)
wilcox.test(dust, dust[-na.ind]) # p-value = 5.991e-05
mis = subset(data, select = c("co2","dust"))
mis$co2[na.ind] = NA
mis$Missing = ifelse(is.na(mis$co2), "Missing", "Observed")
p1 = ggplot(data = mis) + geom_histogram(aes(x=dust,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 20) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of dust (15% missing)") # missing 전후 density 비교

set.seed(100)
na.ind = sample(which(m_sss>quantile(m_sss, 0.5)), nrow(data)*0.3)
wilcox.test(dust, dust[-na.ind]) # p-value < 2.2e-16
mis$co2[na.ind] = NA
mis$Missing = ifelse(is.na(mis$co2), "Missing", "Observed")
p2 = ggplot(data = mis) + geom_histogram(aes(x=dust,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 20) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of dust (30% missing)")

set.seed(100)
na.ind = sample(which(m_sss>quantile(m_sss, 0.5)), nrow(data)*0.45)
wilcox.test(dust, dust[-na.ind]) # p-value< 2.2e-16
mis$co2[na.ind] = NA
mis$Missing = ifelse(is.na(mis$co2), "Missing", "Observed")
p3 = ggplot(data = mis) + geom_histogram(aes(x=dust,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 20) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of dust (45% missing)")

require(gridExtra)
grid.arrange(p1, p2, p3, ncol=3) # dust


#humidity missing
shapiro.test(humidity)
cor(humidity, floor, method = 'spearman') # 0.2635021
m_sss = 0.264*floor
m_sss = (m_sss - min(m_sss))/(max(m_sss)-min(m_sss))
set.seed(100)
na.ind = sample(which(m_sss>=quantile(m_sss, 0.5)), nrow(data)*0.45)
p1 =hist(humidity, breaks = 30)
p2 = hist(humidity[-na.ind], breaks = 30)
plot(p1,col=rgb(0,0,1,1/4), main="Histogram of Humidity (30% missing)")
plot(p2, col=rgb(1,0,0,1/4),add=T)

wilcox.test(floor, floor[-na.ind]) # p-value = 2.789e-05
wilcox.test(humidity,humidity[-na.ind])
mis = subset(data, select = c("humidity","floor"))
mis$humidity[na.ind] = NA
mis$Missing = ifelse(is.na(mis$humidity), "Missing", "Observed")
p1 = ggplot(data = mis) + geom_histogram(aes(x=floor,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 5) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of floor (5% missing)") # missing 전후 density 비교


set.seed(100)
na.ind = sample(which(m_sss>=quantile(m_sss, 0.5)), nrow(data)*0.15)

p1 = hist(humidity)
p2 = hist(humidity[-na.ind])
plot(p1,col=rgb(0,0,1,1/4), main="Histogram of Humidity (15% missing)")
plot(p2, col=rgb(1,0,0,1/4),add=T)


wilcox.test(floor,floor[-na.ind]) # p-value = 2.789e-05
wilcox.test(humidity,humidity[-na.ind])
mis$humidity[na.ind] = NA
mis$Missing = ifelse(is.na(mis$humidity), "Missing", "Observed")
p2 = ggplot(data = mis) + geom_histogram(aes(x=floor,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 5) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of floor (30% missing)") # missing 전후 density 비교

set.seed(100)
na.ind = sample(which(m_sss>=quantile(m_sss, 0.5)), nrow(data)*0.45)

p1 = hist(humidity)
p2 = hist(humidity[-na.ind])
plot(p1,col=rgb(0,0,1,1/4), main="Histogram of Humidity (30% missing)")
plot(p2, col=rgb(1,0,0,1/4),add=T)

wilcox.test(floor,floor[-na.ind]) # p-value = 2.789e-05
wilcox.test(humidity,humidity[-na.ind])
mis$humidity[na.ind] = NA
mis$Missing = ifelse(is.na(mis$humidity), "Missing", "Observed")
p3 = ggplot(data = mis) + geom_histogram(aes(x=floor,fill = Missing, y = ..density..), alpha = 0.2, color = 1, bins = 5) +
  scale_fill_manual(values=c("blue","red")) +
  ggtitle("Density plot of floor (45% missing)") # missing 전후 density 비교
require(gridExtra)
grid.arrange(p1, p2, p3, ncol=3)



#######################################################################################
#15% missing

make_NA=function(df=data,var='SSS', hum=0,temp=0.376, dust_coef=0, f_coef=0, perc=0.15){
  m_sss = hum*humidity+temp*temperature+dust_coef*dust+f_coef*floor
  m_sss = (m_sss - min(m_sss))/(max(m_sss)-min(m_sss))
  set.seed(100)
  na.ind = sample(which(m_sss>=quantile(m_sss, 0.5)), nrow(data)*perc)
  print(length(na.ind))
  df[,var][na.ind]=NA
  return(df)
}
real_final=function(perc1=0.15, perc2=0.05){
  mis_SSS=make_NA(df=data, perc=perc1)
  mis_floor=make_NA(mis_SSS, 'humidity',0,0,0,0.264,perc2)
  mis_co2=make_NA(mis_floor, 'co2',0,0,0.341,0,perc1 )
  final=mis_fire=make_NA(mis_co2, 'fire',-0.374,0.453,0,0,perc1)
  return(final)
}

#15% missing
data15_h5=real_final(0.15, 0.05)
data15_h15=real_final(0.15, 0.15)
data15_h45=real_final(0.15, 0.45)

cor_order=function(data=data){
  df = sapply(data.frame(cor(data[complete.cases(data),])),abs)
  df2 = sapply(data.frame(cor(data[complete.cases(data),], method="spearman")),abs)
  order = colMeans(df)
  order[8] = colMeans(df2)[8]
  return(order)
}
cor_order(data15_h45)
d1 = data[complete.cases(data15_h45),]
wilcox.test(humidity, d1$humidity)


length(d1$humidity)

#30% missing
data30_h5=real_final(0.3, 0.05)
data30_h15=real_final(0.3, 0.15)
data30_h45=real_final(0.3, 0.45)
colSums(is.na(data15_h5))

d1 = data[complete.cases(data30_h45),]
wilcox.test(humidity, d1$humidity)

#45% missing
data45_h5=real_final(0.45, 0.05)
data45_h15=real_final(0.45, 0.15)
data45_h45=real_final(0.45, 0.45)
colSums(is.na(data45_h45))

cor_order(data45_h45)

d1 = data[complete.cases(data45_h5),]
p1 = hist(humidity,  breaks = 25)
p2 = hist(d1$humidity, breaks = 25)
plot(p1,col="purple", main="Histogram of 5% missing humidity(30%)")
plot(p2, col="pink",add=T)


p1 = hist(humidity,  breaks = 25)
p2 = 
plot(p1,col="purple", main="Histogram of 5% missing humidity(30%)")
points(plot(p2, col="pink",add=T))


library(mice)
library(VIM)
aggr(data30_h45, col=mdc(1:2), numbers=TRUE, sortVars=TRUE, labels=names(data15_h5), cex.axis=.7, gap=3, ylab=c("Proportion of missingness","Missingness Pattern of 45%"))
