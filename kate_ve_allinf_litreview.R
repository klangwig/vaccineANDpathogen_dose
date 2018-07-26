#kate_ve_allinf_litreview

rm(list=ls()) # clears workspace

###load
library(ggplot2)
require(deSolve)
library(lme4)
library(plotrix)
library(ggplot2)
library(reshape2)
library(reshape)
library(brglm)
library(MASS)
library(boot)
library(bbmle)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(tidyverse)

#read in data
sysdat = read_csv("sysdat.csv")

sysdat.inf = sysdat %>%
  group_by(disease)%>%
  filter(response_variable=="infection")

sysdat.sum = sysdat %>%
  group_by(disease,group)%>%
  summarise(length(unique(dose)))
View(sysdat.sum)

sysdat.inf.sum = sysdat.inf %>%
  group_by(disease,group)%>%
  summarise(length(unique(dose)))
unique(sysdat.inf.sum$disease)

#i used to filter to diseases with at least 3 doses, but look at any disease!

##use these disease for illustrating effect
unique(sysdat.inf$disease)

##summary figures of prevalence by disease type and pathogen dose
sysdat.inf$frac=sysdat.inf$pos/sysdat.inf$tot

r=ggplot(sysdat.inf, aes(x=log10(dose*0.001),y=frac,col=group,shape=as.factor(vax.dose)))+#,shape=as.factor(vax.dose)
  geom_point(size=2)+
  geom_line()+
  facet_wrap(~disease, scales = "free_x")+
  #geom_line(data=data4,aes(x=log10(dose), y=Prevalence, col=Group),size=1)+#,col="red"
  ylab("Fraction Infection/Clinical")+
  xlab(expression(log[10]~Exposure~Concentration))+
  #  coord_cartesian(ylim=c(-0.1,.65))
  theme_bw() +
  theme(axis.title=element_text(size=20),axis.text.y=element_text(size=15),axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.title=element_blank(),legend.text = element_text(size=20,face="italic"))
r

#literature review - using diseases with >3 dose groups and infection as outcome variable
sysdat.inf$trmt.type = sysdat.inf$group
unique(sysdat.inf$disease)

##malaria#
sysdat.inf[sysdat.inf$disease=="malaria",]
malaria = sysdat.inf[sysdat.inf$disease=="malaria",]
malaria$bn = cut(malaria$dose,breaks = c(0,1,2,4,8,12,20,32),include.lowest = T,labels = F)
pos = aggregate(pos~group+disease+bn+organism+study+response_variable+trmt.type,FUN=sum,data=malaria)
Ns = aggregate(tot~group+disease+bn+organism+study+response_variable+trmt.type,FUN=sum,data=malaria)
pos$tot = Ns$tot
pos$frac = pos$pos / pos$tot
pos$neg = pos$tot - pos$pos
colnames(pos)[3]="dose"
pos
#take out old malaria
sysdat.inf = subset(sysdat.inf, disease!="malaria")
#bring in new malaria
sysdat.inf = bind_rows(sysdat.inf,pos)
#malaria is good
unique(sysdat.inf$disease)

##mareks#
sysdat.inf[sysdat.inf$disease=="mareks",]
mareks = sysdat.inf[sysdat.inf$disease=="mareks",]
mareks2=aggregate(pos~group+disease+dose+organism+study+response_variable+trmt.type,FUN=sum,data=mareks);head(mareks2)
marekscount=aggregate(tot~group+disease+dose+organism+study+response_variable+trmt.type,FUN=sum,data=mareks);head(marekscount)
mareks2$tot=marekscount$tot
mareks2$neg=mareks2$tot-mareks2$pos
mareks2$frac=mareks2$pos/mareks2$tot
mareks2
#take out old mareks
sysdat.inf = subset(sysdat.inf, disease!="mareks")
#bring in new mareks
sysdat.inf = bind_rows(sysdat.inf,mareks2)
unique(sysdat.inf$disease)

#salmonella
sysdat.inf[sysdat.inf$disease=="Salmonella typhimurium",]
typh = sysdat.inf[sysdat.inf$disease=="Salmonella typhimurium",]
typh2=aggregate(pos~group+disease+dose+organism+study+response_variable+trmt.type,FUN=sum,data=typh);head(typh2)
typhcount=aggregate(tot~group+disease+dose+organism+study+response_variable+trmt.type,FUN=sum,data=typh);head(typhcount)
typh2$tot=typhcount$tot
typh2$neg=typh2$tot-typh2$pos
typh2$frac=typh2$pos/typh2$tot
#take out old typh
sysdat.inf = subset(sysdat.inf, disease!="Salmonella typhimurium")
#bring in new typh
sysdat.inf = bind_rows(sysdat.inf,typh2)
unique(sysdat.inf$disease)


##brucella##
sysdat.inf[sysdat.inf$disease=="Brucella abortus",]
unique(sysdat.inf$group[sysdat.inf$disease=="Brucella abortus"])
#brucella good

#polio#
polio = sysdat.inf[sysdat.inf$disease=="polio",]
polio2=polio[polio$group=="B"|polio$group=="D",]
polio2$trmt.type[polio2$group=="B"]="control"
polio2$trmt.type[polio2$group=="D"]="vaccine"
#take out old polio
sysdat.inf = subset(sysdat.inf, disease!="polio")
#bring in new polio
sysdat.inf = bind_rows(sysdat.inf,polio2)

##reovirus##
unique(sysdat.inf$group[sysdat.inf$disease=="reovirus"])
#looks good

#tick-borne encephalitis
unique(sysdat.inf$group[sysdat.inf$disease=="tick-borne encephalitis"])
#tbe is okay although remember this is clinical signs

unique(sysdat.inf$trmt.type)
sysdat.inf = sysdat.inf[sysdat.inf$trmt.type=="control"|sysdat.inf$trmt.type=="vaccine",]
#View(sysdat.inf)

####nice combined figure of 9 diseases#####
xdat=as.data.frame(matrix(nrow=1,ncol=1))
sysdat.inf$frac=sysdat.inf$pos/sysdat.inf$tot

test1=dcast(sysdat.inf,dose+disease~trmt.type, value.var="frac",drop=T, fun.aggregate = mean)
names(test1)=c("dose","disease","control.frac","vaccine.frac")
test2=dcast(sysdat.inf,dose+disease~trmt.type, value.var="pos",drop=T,fun.aggregate = sum)
names(test2)=c("dose","disease","control.pos","vaccine.pos")
test3=dcast(sysdat.inf,dose+disease~trmt.type, value.var="tot",drop=T, fun.aggregate = sum)
names(test3)=c("dose","disease","control.tot","vaccine.tot")
sysdat2=test1
sysdat2$control.pos=test2$control.pos;sysdat2$vaccine.pos=test2$vaccine.pos;
sysdat2$control.tot=test3$control.tot;sysdat2$vaccine.tot=test3$vaccine.tot;
##sysdat2- this is the wide format version which is easier to work with for VE - by brute force...
head(sysdat2)

sysdat2$rr=sysdat2$vaccine.frac/sysdat2$control.frac
sysdat2$ve=1-sysdat2$rr
sysdat2$rr2=sysdat2$rr
sysdat2$rr2[sysdat2$ve==1]=0.1
sysdat2$vaccine.pos2=sysdat2$vaccine.pos
sysdat2$vaccine.pos2[sysdat2$vaccine.pos==0]=1
sysdat2$errs=1/(sysdat2$vaccine.pos2)+1/sysdat2$control.pos-1/sysdat2$vaccine.tot-1/sysdat2$control.tot
sysdat2$se=sqrt(sysdat2$errs)

sysdat2$lrr=log(sysdat2$rr2)
sysdat2$Lrrci=exp(sysdat2$lrr-1.96*sysdat2$se)
sysdat2$Urrci=exp(sysdat2$lrr+1.96*sysdat2$se)

#polio$VEese=1-exp(polio$se)
#with 95% confidence, the rel. risk of getting ihnv w/ vaccine is LCI to UCI of the control
sysdat2$Lveci=1-sysdat2$Lrrci
sysdat2$Uveci=1-sysdat2$Urrci
sysdat2$Lveci[sysdat2$ve==1]=1

lifeisbeautiful=c("#801637","#047878","#FFB733","#F57336","#C22121")
flatdesign=c("#334D5C","#45B29D","#EFC94C","#E27A3F","#DF4949")
sysdat2$disease2 = NA
sysdat2$disease2=as.character(sysdat2$disease)
sysdat2$disease2[sysdat2$disease=="polio"]="Poliovirus"
sysdat2$disease2[sysdat2$disease=="malaria"]="Malaria"
sysdat2$disease2[sysdat2$disease=="tick-borne encephalitis"]="Tick-borne encephalitis"
sysdat2$disease2[sysdat2$disease=="reovirus"]="Reovirus"
sysdat2$disease2[sysdat2$disease=="dhbv"]="Duck Hepatitis B"
sysdat2$disease2[sysdat2$disease=="mareks"]="Mareks Disease"

unique(sysdat2$disease2)
library(viridis)
r2=ggplot(sysdat2, aes(x=log10(dose),y=ve,col=disease2))+
  geom_point(size=2)+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=Lveci, ymax=Uveci), width=.1,color="gray",linetype="dotted")+
  #scale_color_manual(values=flatdesign)+
  scale_color_viridis(discrete=T)+
  facet_wrap(~disease2, scales="free")+
  ylab("Protection Efficacy +/- 95% CI")+
  xlab(expression(log[10]~Dose))+
  coord_cartesian(ylim=c(-0.01,1.01))+
  theme_bw() +
  theme(axis.title=element_text(size=20),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        strip.background = element_rect(fill="gray97"),
        strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position="none",
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"),
        #legend.direction = "vertical",
        legend.box="horizontal")

r2

ggsave(file="/Users/klangwig/Dropbox/fish vaccines/figs/alldiseases_inf_ve_by_dose_20JUL2018.pdf",width=13,height=6,units="in",limitsize=F,useDingbats=FALSE) #saves g

##positives
sysdat.inf$disease2 = NA
sysdat.inf$disease2=as.character(sysdat.inf$disease)
sysdat.inf$disease2[sysdat.inf$disease=="polio"]="Poliovirus"
sysdat.inf$disease2[sysdat.inf$disease=="malaria"]="Malaria"
sysdat.inf$disease2[sysdat.inf$disease=="tick-borne encephalitis"]="Tick-borne encephalitis"
sysdat.inf$disease2[sysdat.inf$disease=="reovirus"]="Reovirus"
sysdat.inf$disease2[sysdat.inf$disease=="dhbv"]="Duck Hepatitis B"
sysdat.inf$disease2[sysdat.inf$disease=="mareks"]="Mareks Disease"

sysdat.inf$trmt.type2[sysdat.inf$trmt.type=="control"]="Control"
sysdat.inf$trmt.type2[sysdat.inf$trmt.type=="vaccine"]="Vaccine/Antibody"

r=ggplot(sysdat.inf, aes(x=log10(dose),y=frac,col=trmt.type2))+#,shape=as.factor(vax.dose)
  geom_point(size=2)+
  geom_line()+
  facet_wrap(~disease2, scales = "free_x")+
  #geom_line(data=data4,aes(x=log10(dose), y=Prevalence, col=Group),size=1)+#,col="red"
  ylab("Fraction Infected/Clinical")+
  xlab("")+
  #xlab(expression(log[10]~Exposure~Concentration))+
  scale_color_manual(values=c("red","blue"))+
  theme_bw() +
  theme(axis.title=element_text(size=20),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        strip.background = element_rect(fill="gray97"),
        strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position="top",
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"),
        #legend.direction = "vertical",
        legend.box="horizontal")
r

ggsave(file="/Users/klangwig/Dropbox/fish vaccines/figs/alldiseases_infection_by_dose_litreview_20JUL2018.pdf",width=13,height=6,units="in",limitsize=F,useDingbats=FALSE) #saves g

g2 <- arrangeGrob(r, r2, nrow=2) #generates g2

ggsave(file="/Users/klangwig/Dropbox/fish vaccines/figs/combined_litreview.pdf", g2,width=13,height=12,units="in",limitsize=F,useDingbats=FALSE) #saves g


##statistical analyses
head(sysdat2); head(sysdat.inf)
#ordinal regression
#Is it worth formalizing the decline by logistic regression 
#logit Pr(inf|dose) = b0+b1*logdose +b2*vax+b3*logdose*vax and testing significance of b3? 
#It is crude but maybe a good idea? Or maybe just showing the CI on VE so reader can decide if s/he believes it?
sysdat.inf2 = sysdat.inf %>%
  group_by(disease2,trmt.type) %>%
  mutate(no.doses = length(unique(dose)))%>%
  arrange(dose) %>%
  mutate(dose.o = row_number())



m <-glmer(frac ~ dose.o*trmt.type + (1|disease2), weights = tot, data = sysdat.inf2, family="binomial")
summary(m)
library(car)
Anova(m)

head(sysdat2)

sysdat3 = sysdat2 %>%
  group_by(disease2) %>%
  mutate(no.doses = length(unique(dose)))%>%
  arrange(dose) %>%
  mutate(dose.o = row_number())

m2 <-lmer(ve ~ dose.o + (1|disease2), data = sysdat3)
summary(m2)


