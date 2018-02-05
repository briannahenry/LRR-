library(brms)

#trt2=atrazine; trt3=selenium; trt4=combined#

dbentab$trt0f<-as.factor(dbentab$trt)
dbentab$date0f<-as.factor(dbentab$date)


#######Total Larval Abundance Model and Plot##########
get_prior(nobugs~date0f+trt0f+date0f*trt0f +(1|tank),data=dbentab,family=poisson(link="log"))

larvabund<-brm(nobugs~date0f+trt0f+date0f*trt0f +(1|tank),data=dbentab,family=poisson(link="log"),
           prior=c(prior(cauchy(0,1),class=sd),
                   prior(normal(0,10),class=b)))

print(larvabund)
pp_check(larvabund,type="hist") 
pp_check(larvabund,type="boxplot")

marginal_effects(larvabund)

laplot<-marginal_effects(larvabund, robust=FALSE)
laplot<-as.data.frame(laplot$`date0f:trt0f`)
laplot$date2<-factor(laplot$date0f, levels = c("6/20/2016","7/6/2016","8/2/2016"))
laplot$est2<-(laplot$estimate__)*31.25
laplot$upper2<-(laplot$upper__)*31.25
laplot$lower2<-(laplot$lower__)*31.25
laplot$trt2<-factor(laplot$trt0f, levels = c("Control", "Atrazine", "Selenium", "Atrazine and Selenium"))



pd<- position_dodge(0.4) # stagger points .2 to the left and right

larvabundplot<-ggplot(laplot,aes(x=date2,y=est2, color=trt0f))+
  geom_errorbar(aes(ymin=lower2, ymax=upper2), width=.1, position=pd)+
  geom_point(position=pd, size=4)+
  xlab("Date Sampled")+
  ylab("Total Larval Abundance (#/m2)")+
  scale_color_manual(name="Treatment", labels = c("Control", "Atrazine", "Selenium", "Atrazine and Selenium"),  values = c("grey0","grey45","grey65","grey80"))+
  theme_classic()+
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black"))+
  theme(axis.title.x = element_text(face="bold"))+
  theme(axis.title.y = element_text(face="bold"))+
  theme(legend.title = element_text(face="bold", size="11"))

larvabundplot

####Labundbytrt<-ggplot(laplot,aes(x=trt2,y=est2))+
  geom_errorbar(aes(ymin=lower2, ymax=upper2), width=.1)+
  geom_point(size=4)+
  xlab("Date Sampled")+
  ylab("Total Larval Abundance (#/m2)")+
  theme_classic()+
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black"))+
  theme(axis.title.x = element_text(face="bold"))+
  theme(axis.title.y = element_text(face="bold"))+
  theme(legend.title = element_text(face="bold", size="11"))####



ggsave("larvabundplot.tiff", larvabundplot, dpi=800)

############Date/Time Differences in Posterior#######
mlpost<-posterior_samples(larvabund)
str(mlpost)


ltrt1_1<-mlpost$b_Intercept 
ltrt1_2<-mlpost$b_Intercept+mlpost$b_date0f2
ltrt1_3<-mlpost$b_Intercept+mlpost$b_date0f3

ltrt2_1<-mlpost$b_Intercept+mlpost$b_trt0f2 
ltrt2_2<-mlpost$b_Intercept+mlpost$b_date0f2+mlpost$b_trt0f2+mlpost$`b_date0f2:trt0f2`
ltrt2_3<-mlpost$b_Intercept+mlpost$b_date0f3+mlpost$b_trt0f2+mlpost$`b_date0f3:trt0f2`

ltrt3_1<-mlpost$b_Intercept+mlpost$b_trt0f3 
ltrt3_2<-mlpost$b_Intercept+mlpost$b_date0f2+mlpost$b_trt0f3+mlpost$`b_date0f2:trt0f3`
ltrt3_3<-mlpost$b_Intercept+mlpost$b_date0f3+mlpost$b_trt0f3+mlpost$`b_date0f3:trt0f3`

ltrt4_1<-mlpost$b_Intercept+mlpost$b_trt0f4
ltrt4_2<-mlpost$b_Intercept+mlpost$b_date0f2+mlpost$b_trt0f4+mlpost$`b_date0f2:trt0f4`
ltrt4_3<-mlpost$b_Intercept+mlpost$b_date0f3+mlpost$b_trt0f4+mlpost$`b_date0f3:trt0f4`

#Differences on Date 1
#Trt 1 & 2
ltrt1_2_d1<-(exp(ltrt1_1)*31.25)-(exp(ltrt2_1)*31.25)
mean(ltrt1_2_d1)
quantile(ltrt1_2_d1,probs=c(0.025,0.975))
sum(ltrt1_2_d1>0)/4000 #72.3% prob diff > 0

#Trt 1 & 3
ltrt1_3_d1<-(exp(ltrt1_1)*31.25)-(exp(ltrt3_1)*31.25)
mean(ltrt1_3_d1)
quantile(ltrt1_3_d1,probs=c(0.025,0.975))
sum(ltrt1_3_d1>0)/4000 #80.8% prob diff > 0

#Trt 1 & 4
ltrt1_4_d1<-(exp(ltrt1_1)*31.25)-(exp(ltrt4_1)*31.25)
mean(ltrt1_4_d1)
quantile(ltrt1_4_d1,probs=c(0.025,0.975))
sum(ltrt1_4_d1>0)/4000 #84.3% prob diff > 0

#Differences on Date 2
#Trt 1 & 4
ltrt1_4_d2<-(exp(ltrt4_2)*31.25)-(exp(ltrt1_2)*31.25)
mean(ltrt1_4_d2)
quantile(ltrt1_4_d2,probs=c(0.025,0.975))
sum(ltrt1_4_d2>0)/4000 #84.5% prob diff > 0

#Trt 2 & 4 
ltrt2_4_d2<-(exp(ltrt4_2)*31.25)-(exp(ltrt2_2)*31.25)
mean(ltrt2_4_d2)
quantile(ltrt2_4_d2,probs=c(0.025,0.975))
sum(ltrt2_4_d2>0)/4000 #95% prob diff > 0

#Trt 3 & 4 
ltrt3_4_d2<-(exp(ltrt4_2)*31.25)-(exp(ltrt3_2)*31.25)
mean(ltrt3_4_d2)
quantile(ltrt3_4_d2,probs=c(0.025,0.975))
sum(ltrt3_4_d2>0)/4000 #99.5% prob diff > 0

#Trt 1 & 3
ltrt1_3_d2<-(exp(ltrt1_2)*31.25)-(exp(ltrt3_2)*31.25)
mean(ltrt1_3_d2)
quantile(ltrt1_3_d2,probs=c(0.025,0.975))
sum(ltrt1_3_d2>0)/4000 #94.5% prob diff > 0

#Differences on Date 3
#Trt 1 & 2
ltrt1_2_d3<-(exp(ltrt1_3)*31.25)-(exp(ltrt2_3)*31.25)
mean(ltrt1_2_d3)
quantile(ltrt1_2_d3,probs=c(0.025,0.975))
sum(ltrt1_2_d3>0)/4000 #87% prob diff > 0

#Trt 1 & 3
ltrt1_3_d3<-(exp(ltrt1_3)*31.25)-(exp(ltrt3_3)*31.25)
mean(ltrt1_3_d3)
quantile(ltrt1_3_d3,probs=c(0.025,0.975))
sum(ltrt1_3_d3>0)/4000 #98.9% prob diff > 0

#Trt 1 & 4
ltrt1_4_d3<-(exp(ltrt1_3)*31.25)-(exp(ltrt4_3)*31.25)
mean(ltrt1_4_d3)
quantile(ltrt1_4_d3,probs=c(0.025,0.975))
sum(ltrt1_4_d3>0)/4000 #98% prob diff > 0
