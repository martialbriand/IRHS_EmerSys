
#title: "Single seed microbiota: assembly and transmission from parent plant to seedling"
#author: "Chesneau et al., 2021"
#date: "21/05/2021"

setwd("~/These2/Script R/F2S-diversite/F2S-year-4/Markdown/Dossier_run_script_vf/Analyse_beatrice/")
library(stats)
library(data.table)
library(MASS)
library(bvls)
library(BiocManager)
library(TailRank)
library(DescTools)
library(fitdistrplus)
plot.new()

################################
# USER SUPPLIED INFORMATIONS
# abundance file
ficname="otu_table_16Sb.csv"
sep=','
Nplanttot=24
# to simulate the neutral model: we need the number of draws in the distribution of N (number of simulated seeds) and
# a grid for the model parameters S and R
Ndraw=1000
S=2:150
R= (1:250)*0.000001#mu/((1-mu)*Cstar)
#data filtering
type=1 #1:remove individual OTUs whose proportion within one seed is less than threshold (0.3% works, 0.5% works fine)
#2: remove all OTUs whose cumulated proportions within one seed is less than threshold (1% works fine)
threshold=1. #threshold for data filtering in % 
##################################

#read data files, samples are lines, ASV are columns
#select Beans if necessary, list sampling dates
ASVbean<-read.csv(file = ficname,sep=sep,header=TRUE,)
#ASVbean=ASVdata[which(ASVdata$Plant.species=="Bean"),] 
day=unique(ASVbean$DAP)

indexASV=which(substr(colnames(ASVbean),1,3)=="ASV")
C=ASVbean[,indexASV]
# how does the distribution of the number of ASV detected in seeds look like?
#the answer is: almost uniform 2000-23700
par(mfrow=c(1,1))
Ntot=rowSums(C)
plot(rownames(table(rowSums(C)),cumsum(table(rowSums(C)))))
x=runif(2000, min=min(rowSums(C)), max=max(rowSums(C)))
qqplot(x,rowSums(C))
qqline(rowSums(C), distribution = function(p) qunif(p, min=min(rowSums(C)), max=max(rowSums(C))), probs = c(0., 1), col = 2)
#just to know


# now we build the empirical distribution per date and plant
# and fit an uniform on the Ntot per date

dayplant=NULL
Fday=list()
M1list=list()
M2list=list()
M3list=list()
nl=0
lambda=NULL
sd=NULL
res=NULL
Filtering=NULL

#loop on dates
for (j in day) {
  indexDAP=which(ASVbean$DAP==j)
  ASVbeanj=ASVbean[indexDAP,]
  Ntotj=Ntot[indexDAP]
  qu=quantile(Ntotj,probs=c(0.1,0.9))
  lambda=rbind(lambda, qu)
#
#now loop on plants within each day
  plant=unique(ASVbeanj$Plant)
  for (p in plant){
    indexP=which(ASVbeanj$Plant==p)
    ASVbeanjp=ASVbeanj[indexP,]
    nseed=dim(ASVbeanjp)[1]
# number of seeds per plant and date
    A=transpose(ASVbeanjp[,indexASV])
    row.names(A)=colnames(ASVbeanjp[indexASV])
    index=which(rowSums(A)>0)
    Splant=sum(index>0)
    data0=sum(A>0)
# number of ASV present at least once per time and date
    A=A[index,]
    L=sort(A[,1], decreasing = TRUE,index.return=TRUE)
    B=as.vector(L$x)
    if (type==2){
    Nthresh=floor(sum(B)*(1-threshold/100))
    filter=which(cumsum(B)>Nthresh)
    B[filter[2:length(filter)]]=0 }
    else {
    Nthresh=floor(sum(B)*threshold/100)
    filter=which(B<=Nthresh)
    B[filter]=0 }
    #here we need sorting
    F=B/sum(B) #frequency per seed
    for (i in (2:nseed)){
      l=sort(A[,i], decreasing = TRUE, index.return=TRUE)
      b=as.vector(l$x)
      if (type==2){
      Nthresh=floor(sum(b)*(1-threshold/100))
      filter=which(cumsum(b)>Nthresh)
      b[filter[2:length(filter)]]=0}
      else {
        Nthresh=floor(sum(b)*threshold/100)
        filter=which(b<=Nthresh)
        b[filter]=0 }
      B=cbind(B,b)
      F=cbind(F,b/sum(b))
    }
  # now compute and store first three moment on the non zero asv
    m1Fpos=0*(1:dim(F)[2])
    m2Fpos=0*(1:dim(F)[2])
    m3Fpos=0*(1:dim(F)[2])
    for (s in 1:dim(F)[2]){
      m1Fpos[s]=mean(F[which(F[,s]>0),s])
      m2Fpos[s]=mean(F[which(F[,s]>0),s]^2)
      m3Fpos[s]=mean(F[which(F[,s]>0),s]^3)
    }
    datafiltered=sum(B>0)
    Filtering=c(Filtering, datafiltered/data0)
    nl=nl+1
    Fday[[nl]]=F
    M1list[[nl]]=m1Fpos
    M2list[[nl]]=m2Fpos
    M3list[[nl]]=m3Fpos
    daynumber=which(day==j)
    plantnumber=which(plant==p)
    dayplant[[nl]]=c(daynumber,plantnumber)
  }
}
# uncomment to see the impact of filtering: proportion of non zero abundance data kept after filtering for each plant
# hist(Filtering)

# moment estimation of S and R related to neutral model
# we simulate the neutral model and compare with the true distribution
# 3 slightly different weighted criteria to be minimized, see text. basically we use popt
popt=NULL
popt2=NULL
popt3=NULL
#
for (nl in 1:Nplanttot){
p=dayplant[[nl]]
set.seed(nl)
Nsample=sample(lambda[p[1],1]:lambda[p[1],2], Ndraw, replace=T)
# draw values for N in the suitable uniform distribution 
err=NULL
err2=NULL
err3=NULL
param=NULL 
count=0
M1obs=mean(M1list[[nl]])
M2obs=mean(M2list[[nl]])
M3obs=mean(M3list[[nl]])
for (s in S) {
   for (r in R) {
count=count+1
Ic=r*(Nsample-1)
alphac=Ic/s
betac=Ic-alphac
#compute parameters for the beta-binomial
#
P0=dbb(0,Nsample,alphac,betac)
# probability of zero to threshold under the beta-binomial model
Em1=mean(1/(s*(1-P0)))
# marginal (on N) first moment of bb distribution conditional on n>0
errM1=abs(Em1-M1obs)
numerator=1+alphac+betac/Nsample
denominator1=(1+Ic)*s*(1-P0)
Em2=mean(numerator/denominator1)
# marginal (on N) second moment of bb distribution conditional on n>0
errM2=abs(Em2-M2obs)
numerator=(1+alphac)*((2+alphac)+3*betac/Nsample)+betac*(betac-alphac)/(Nsample^2)
denominator2=denominator1*(2+Ic)
Em3=mean(numerator/denominator2)
# marginal (on N) third moment of bb distribution conditional on n>0
errM3=abs(Em3-M3obs)
param[[count]]=c(s,r)
err3=c(err3, sqrt(errM1^2+(errM2/M1obs)^2+(errM3/M2obs)^2))
err2=c(err2, sqrt((errM1/M1obs)^2+(errM2/M2obs)^2+(errM3/M3obs)^2))
err=c(err, sqrt(errM1^2+(errM2)^2+(errM3)^2))
# values for the three criteria
      }
   }
bestp=param[[which(err==min(err))]]
ic=bestp[2]*(mean(Nsample)-1) #R*(expected (N)-1)
alpha=ic/bestp[1]# parameter of the bb distribution
popt=rbind(popt,c(p,param[[which(err==min(err))]], ic, alpha, err[which(err==min(err))]))
bestp=param[[which(err2==min(err2))]]
ic=bestp[2]*(mean(Nsample)-1) #R*(expected (N)-1)
alpha=ic/bestp[1]
popt2=rbind(popt2,c(p,param[[which(err2==min(err2))]], ic, alpha, err2[which(err2==min(err2))]))
bestp=param[[which(err3==min(err3))]]
ic=bestp[2]*(mean(Nsample)-1) #R*(expected (N)-1)
alpha=ic/bestp[1]
popt3=rbind(popt3,c(p,param[[which(err3==min(err3))]], ic, alpha, err3[which(err3==min(err3))]))
}

#some plots
par(mfrow=c(1,2))
boxplot(popt[,7]~popt[,1],outline=FALSE,names = day,main = "Error")
boxplot(popt[,3]~popt[,1],outline=FALSE,names = day,main = "S")
par(mfrow=c(1,2))
boxplot(popt[,5]~popt[,1],outline=FALSE,names = day,main = "ratio Iplant")
boxplot(popt[,4]~popt[,1],outline=FALSE,names = day,main = "R")


#Write file
write.csv(popt, "popt_b_16S")

#Prepare data plot
poptf <- data.frame(popt)
colnames(poptf) <- c("DAP", "X2", "S", "R", "ratio_Iplant", "X6", "Error")
poptf$DAP<-ifelse(poptf$DAP == 1,'D27',
                  ifelse(poptf$DAP == 2,'D33',
                         ifelse(poptf$DAP == 3,'D42',
                                ifelse(poptf$DAP == 4,'D50','D50'))))

#Prepare palette color for habitats
name_DAP <- c("D02","D27", "D33", "D42", "D50", "D03","D26", "D37", "D65")

# desired color palette
DAPPalette <- c( "#F0FEEF","#F5F6CE", "#D0F5A9", "#5FB404", "#0B610B","#FFF8EC", "#F6E3CE", "#FAAC58", "#DF3A01")

# associated EnvType level
names(DAPPalette) <- name_DAP

#Write file
write.csv(poptf, "poptf_b_16S")

#Open data
poptf_b_16S <- read.csv("poptf_b_16S", sep = ",",  check.names=FALSE, row.names=1)
poptf_b_gyrB <- read.csv("poptf_b_gyrB", sep = ",",  check.names=FALSE, row.names=1)

poptf_b_16S$Plant.species <- c("Bean")
poptf_b_gyrB$Plant.species <- c("Bean")

#Plot

library(ggplot2)
theme_set(theme_bw())

Plot_16Sb_Iplant <- ggplot(data=poptf_b_16S, aes_string(x='DAP',y='ratio_Iplant', fill='DAP')) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=DAP)) +
  theme(strip.background = element_rect(color="white", fill="#585858", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "white"),
        axis.title.x = element_text( size = 12, face = "bold"),
        axis.title.y = element_text( size = 12, face = "bold"),
        plot.title = element_text(size = 16),
        axis.text.x = element_text(), 
        axis.text.y = element_text())+ 
  ylab ("Ratio I plant") + xlab ("Day after polination") +
  facet_grid(cols = vars(Plant.species), scales = "free") +
  ggtitle("Ratio I plant - 16S") + 
  scale_fill_manual(values=DAPPalette) + 
  scale_color_manual(values=DAPPalette)
Plot_16Sb_Iplant

Plot_16Sb_R <- ggplot(data=poptf_b_16S, aes_string(x='DAP',y='R', fill='DAP')) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=DAP)) +
  theme(strip.background = element_rect(color="white", fill="#585858", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "white"),
        axis.title.x = element_text( size = 12, face = "bold"),
        axis.title.y = element_text( size = 12, face = "bold"),
        plot.title = element_text(size = 16),
        axis.text.x = element_text(), 
        axis.text.y = element_text())+ 
  ylab ("R") + xlab ("Day after polination") +
  facet_grid(cols = vars(Plant.species), scales = "free") +
  ggtitle("R - 16S") + 
  scale_fill_manual(values=DAPPalette) + 
  scale_color_manual(values=DAPPalette)
Plot_16Sb_R

Plot_16Sb_Error <- ggplot(data=poptf_b_16S, aes_string(x='DAP',y='Error', fill='DAP')) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=DAP)) +
  theme(strip.background = element_rect(color="white", fill="#585858", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "white"),
        axis.title.x = element_text( size = 12, face = "bold"),
        axis.title.y = element_text( size = 12, face = "bold"),
        plot.title = element_text(size = 16),
        axis.text.x = element_text(), 
        axis.text.y = element_text())+ 
  ylab ("Error") + xlab ("Day after polination") +
  facet_grid(cols = vars(Plant.species), scales = "free") +
  ggtitle("Error - 16S") + 
  scale_fill_manual(values=DAPPalette) + 
  scale_color_manual(values=DAPPalette)
Plot_16Sb_Error

Plot_16Sb_S <- ggplot(data=poptf_b_16S, aes_string(x='DAP',y='S', fill='DAP')) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=DAP)) +
  theme(strip.background = element_rect(color="white", fill="#585858", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "white"),
        axis.title.x = element_text( size = 12, face = "bold"),
        axis.title.y = element_text( size = 12, face = "bold"),
        plot.title = element_text(size = 16),
        axis.text.x = element_text(), 
        axis.text.y = element_text())+ 
  ylab ("S") + xlab ("Day after polination") +
  facet_grid(cols = vars(Plant.species), scales = "free") +
  ggtitle("S - 16S") + 
  scale_fill_manual(values=DAPPalette) + 
  scale_color_manual(values=DAPPalette)
Plot_16Sb_S


library(cowplot); packageVersion('cowplot')
Plot_16S_b <- plot_grid(Plot_16Sb_Iplant, Plot_16Sb_R, Plot_16Sb_Error, Plot_16Sb_S)



Plot_gyrBb_Iplant <- ggplot(data=poptf_b_gyrB, aes_string(x='DAP',y='ratio_Iplant', fill='DAP')) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=DAP)) +
  theme(strip.background = element_rect(color="white", fill="#585858", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "white"),
        axis.title.x = element_text( size = 12, face = "bold"),
        axis.title.y = element_text( size = 12, face = "bold"),
        plot.title = element_text(size = 16),
        axis.text.x = element_text(), 
        axis.text.y = element_text())+ 
  ylab ("Ratio I plant") + xlab ("Day after polination") +
  facet_grid(cols = vars(Plant.species), scales = "free") +
  ggtitle("Ratio I plant - gyrB") + 
  scale_fill_manual(values=DAPPalette) + 
  scale_color_manual(values=DAPPalette)
Plot_gyrBb_Iplant

Plot_gyrBb_R <- ggplot(data=poptf_b_gyrB, aes_string(x='DAP',y='R', fill='DAP')) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=DAP)) +
  theme(strip.background = element_rect(color="white", fill="#585858", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "white"),
        axis.title.x = element_text( size = 12, face = "bold"),
        axis.title.y = element_text( size = 12, face = "bold"),
        plot.title = element_text(size = 16),
        axis.text.x = element_text(), 
        axis.text.y = element_text())+ 
  ylab ("R") + xlab ("Day after polination") +
  facet_grid(cols = vars(Plant.species), scales = "free") +
  ggtitle("R - gyrB") + 
  scale_fill_manual(values=DAPPalette) + 
  scale_color_manual(values=DAPPalette)
Plot_gyrBb_R

Plot_gyrBb_Error <- ggplot(data=poptf_b_gyrB, aes_string(x='DAP',y='Error', fill='DAP')) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=DAP)) +
  theme(strip.background = element_rect(color="white", fill="#585858", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "white"),
        axis.title.x = element_text( size = 12, face = "bold"),
        axis.title.y = element_text( size = 12, face = "bold"),
        plot.title = element_text(size = 16),
        axis.text.x = element_text(), 
        axis.text.y = element_text())+ 
  ylab ("Error") + xlab ("Day after polination") +
  facet_grid(cols = vars(Plant.species), scales = "free") +
  ggtitle("Error - gyrB") + 
  scale_fill_manual(values=DAPPalette) + 
  scale_color_manual(values=DAPPalette)
Plot_gyrBb_Error

poptf_b_gyrB_S <- subset(poptf_b_gyrB, !S %in% c(57))

Plot_gyrBb_S <- ggplot(data=poptf_b_gyrB_S, aes_string(x='DAP',y='S', fill='DAP')) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(color=DAP)) +
  theme(strip.background = element_rect(color="white", fill="#585858", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "white"),
        axis.title.x = element_text( size = 12, face = "bold"),
        axis.title.y = element_text( size = 12, face = "bold"),
        plot.title = element_text(size = 16),
        axis.text.x = element_text(), 
        axis.text.y = element_text())+ 
  ylab ("S") + xlab ("Day after polination") +
  facet_grid(cols = vars(Plant.species), scales = "free") +
  ggtitle("S - gyrB") + 
  scale_fill_manual(values=DAPPalette) + 
  scale_color_manual(values=DAPPalette)
Plot_gyrBb_S


library(cowplot); packageVersion('cowplot')
Plot_gyrB_b <- plot_grid(Plot_gyrBb_Iplant, Plot_gyrBb_R, Plot_gyrBb_Error, Plot_gyrBb_S)





#tests on I
wilcox.test(poptf_b_gyrB[2:12,5],poptf_b_gyrB[13:24,5],alternative="less")
wilcox.test(poptf_b_gyrB[7:12,5],poptf_b_gyrB[13:24,5],alternative="less")
wilcox.test(poptf_b_gyrB[7:12,5],poptf_b_gyrB[13:17,5],alternative="less")
wilcox.test(poptf_b_gyrB[7:12,5],poptf_b_gyrB[18:24,5],alternative="less")
wilcox.test(poptf_b_gyrB[1:6,5],poptf_b_gyrB[7:12,5],alternative="greater")

wilcox.test(poptf_b_gyrB[1:6,3],poptf_b_gyrB[7:12,3],alternative="less")
wilcox.test(poptf_b_gyrB[7:12,3],poptf_b_gyrB[13:17,3],alternative="less")
wilcox.test(poptf_b_gyrB[7:12,3],poptf_b_gyrB[18:24,3],alternative="less")

#tests on I
wilcox.test(poptf_b_16S[2:12,5],poptf_b_16S[13:24,5],alternative="less")
wilcox.test(poptf_b_16S[7:12,5],poptf_b_16S[13:24,5],alternative="less")
wilcox.test(poptf_b_16S[7:12,5],poptf_b_16S[13:17,5],alternative="less")
wilcox.test(poptf_b_16S[7:12,5],poptf_b_16S[18:24,5],alternative="less")
wilcox.test(poptf_b_16S[1:6,5],poptf_b_16S[7:12,5],alternative="greater")

wilcox.test(poptf_b_16S[18:24,7],poptf_b_16S[1:17,7],alternative="less")


pairwise.wilcox.test(poptf_b_gyrB[,5], poptf_b_gyrB[,1], p.adjust.method="BH")
# raté, pour les tests multiples ça ne marche pas

pairwise.wilcox.test(poptf_b_gyrB[,3], poptf_b_gyrB[,1], p.adjust.method="BH")
#mais pour la taille du pool local ça fonctionne presque
#on poole les deux classes du milieu
pooltest=poptf_b_gyrB[,1]
pooltest[7:17]=2
pooltest[18:24]=3
pairwise.wilcox.test(poptf_b_gyrB[,3],pooltest, p.adjust.method="BH")
wilcox.test(poptf_b_gyrB[1:17,3],poptf_b_gyrB[18:24,3],alternative="less")


pairwise.wilcox.test(poptf_b_gyrB[,7], poptf_b_gyrB[,1], p.adjust.method="BH")
wilcox.test(poptf_b_gyrB[1:6,7],poptf_b_gyrB[7:24,7],alternative="greater")
wilcox.test(poptf_b_gyrB[18:24,7],poptf_b_gyrB[1:17,7],alternative="less")

par(mfrow=c(2,4))
# so now we do it for the actual marginal (on N) frequency distribution
#with nice titles
Ndraw=500
for (nl in 1:Nplanttot){
  PD=as.vector(Fday[[nl]])
  titre=paste("Day",day[popt[nl,1]],"Plant",popt[nl,2],sep=" ")
  hist(PD[which(PD>0)],(0:100)/100,freq=FALSE,main=titre,xlab="proportion",ylab="Probability density")
  S=popt[nl,3]
  R=popt[nl,4]
  set.seed(nl)
  Nsample=sample(lambda[popt[nl,1],1]:lambda[popt[nl,1],2], Ndraw, replace=T)
  Nmax=max(Nsample)
  Fn=NULL
  Pn=NULL
  for (N in Nsample) {
    IN=R*(N-1)
    alphaN=IN/S
    betaN=IN-alphaN
    #
    P0=dbb(0,N,alphaN,betaN)
    Pn=append(Pn,dbb(1:N,N,alphaN,betaN)/(1-P0),after=length(Pn))
    Fn=append(Fn,(1:N)/N,after=length(Fn))
  }
  L=sort(Fn,index.return=TRUE)
  Pntrie=Pn[L$ix]
  Dens=0*(1:100)
  for (i in 1:100) {Dens[i]=sum(Pntrie[which(L$x>(i-1)/100 & L$x<=i/100)])/Ndraw}
  lines((1:100)/100, 100*Dens,col="red")
}
par(mfrow=c(1,1))
