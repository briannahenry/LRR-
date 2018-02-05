library(brms)

#trt2=atrazine; trt3=selenium; trt4=combined#
#####EMERGENCE MODEL USING BRMS package####

###create factors as needed
biodivnew$trt0f<-as.factor(biodivnew$trt)
biodivnew$date0f<-as.factor(biodivnew$date)

#######Model and Plot#######

get_prior(totno~date0f+trt0f+date0f*trt0f +(1|tank),data=biodivnew,family=poisson(link="log"))

abundadult<-brm(totno~date0f+trt0f+date0f*trt0f +(1|tank),data=biodivnew,family=poisson(link="log"),
           prior=c(prior(cauchy(0,1),class=sd),
                   prior(normal(0,5),class=b)))

print(abundadult)
pp_check(abundadult,type="hist") 
pp_check(abundadult,type="boxplot")

marginal_effects(abundadult, robust=FALSE)

aaplot<-marginal_effects(abundadult, robust=FALSE)
aaplot<-as.data.frame(aaplot$`date0f:trt0f`)
aaplot$date2<-factor(aaplot$date0f, levels = c("6/20/2016","6/28/2016", "7/4/2016", "7/12/2016", "7/18/2016", "7/26/2016", "8/2/2016"))
aaplot$est2<-(aaplot$estimate__)*2.78/4
aaplot$upper2<-(aaplot$upper__)*2.78/4
aaplot$lower2<-(aaplot$lower__)*2.78/4

pd<- position_dodge(0.4) # stagger points .2 to the left and right

adultabundplot<-ggplot(aaplot,aes(x=date2,y=est2, color=trt0f))+
  geom_errorbar(aes(ymin=lower2, ymax=upper2), width=.1, position=pd)+
  geom_point(position=pd, size=4)+
  xlab("Date Sampled")+
  ylab("Total Adult Abundance (#/m2/day)")+
  scale_color_manual(name="Treatment", labels = c("Control", "Atrazine", "Selenium", "Atrazine and Selenium"),  values = c("grey0","grey45","grey65","grey80"))+
  theme_classic()+
  theme(axis.text.x = element_text(color="black"))+
  theme(axis.text.y = element_text(color="black"))+
  theme(axis.title.x = element_text(face="bold"))+
  theme(axis.title.y = element_text(face="bold"))+
  theme(legend.title = element_text(face="bold", size="11"))

adultabundplot

ggsave("adultabundplot.tiff", adultabundplot, dpi=800)


######Date_trt Comparisons########

#extract posteriors
m2post<-posterior_samples(m2brm)
str(m2post)

#estimate posteriors for each date*trt combination. 

trt1_2<-m2post$b_Intercept+m2post$b_date0f2 #posterior treatment 1 on date 2
trt1_3<-m2post$b_Intercept+m2post$b_date0f3 #posterior treatment 1 on date 3
trt1_4<-m2post$b_Intercept+m2post$b_date0f4 #posterior treatment 1 on date 3
trt1_5<-m2post$b_Intercept+m2post$b_date0f5 #posterior treatment 1 on date 3
trt1_6<-m2post$b_Intercept+m2post$b_date0f6 #posterior treatment 1 on date 3
trt1_7<-m2post$b_Intercept+m2post$b_date0f7 #posterior treatment 1 on date 3

trt2_2<-m2post$b_Intercept+m2post$b_date0f2+m2post$b_trt0f2+m2post$`b_date0f2:trt0f2`
trt2_3<-m2post$b_Intercept+m2post$b_date0f3+m2post$b_trt0f2+m2post$`b_date0f3:trt0f2`
trt2_4<-m2post$b_Intercept+m2post$b_date0f4+m2post$b_trt0f2+m2post$`b_date0f4:trt0f2`
trt2_5<-m2post$b_Intercept+m2post$b_date0f5+m2post$b_trt0f2+m2post$`b_date0f5:trt0f2`
trt2_6<-m2post$b_Intercept+m2post$b_date0f6+m2post$b_trt0f2+m2post$`b_date0f6:trt0f2`
trt2_7<-m2post$b_Intercept+m2post$b_date0f7+m2post$b_trt0f2+m2post$`b_date0f7:trt0f2`

trt3_2<-m2post$b_Intercept+m2post$b_date0f2+m2post$b_trt0f3+m2post$`b_date0f2:trt0f3`
trt3_3<-m2post$b_Intercept+m2post$b_date0f3+m2post$b_trt0f3+m2post$`b_date0f3:trt0f3`
trt3_4<-m2post$b_Intercept+m2post$b_date0f4+m2post$b_trt0f3+m2post$`b_date0f4:trt0f3`
trt3_5<-m2post$b_Intercept+m2post$b_date0f5+m2post$b_trt0f3+m2post$`b_date0f5:trt0f3`
trt3_6<-m2post$b_Intercept+m2post$b_date0f6+m2post$b_trt0f3+m2post$`b_date0f6:trt0f3`
trt3_7<-m2post$b_Intercept+m2post$b_date0f7+m2post$b_trt0f3+m2post$`b_date0f7:trt0f3`

trt4_2<-m2post$b_Intercept+m2post$b_date0f2+m2post$b_trt0f4+m2post$`b_date0f2:trt0f4` #posterior treatment 4 on date 2
trt4_3<-m2post$b_Intercept+m2post$b_date0f3+m2post$b_trt0f4+m2post$`b_date0f3:trt0f4` #posterior treatment 4 on date 3
trt4_4<-m2post$b_Intercept+m2post$b_date0f4+m2post$b_trt0f4+m2post$`b_date0f4:trt0f4` #posterior treatment 4 on date 4
trt4_5<-m2post$b_Intercept+m2post$b_date0f5+m2post$b_trt0f4+m2post$`b_date0f5:trt0f4` #posterior treatment 4 on date 5
trt4_6<-m2post$b_Intercept+m2post$b_date0f6+m2post$b_trt0f4+m2post$`b_date0f6:trt0f4` #posterior treatment 4 on date 6
trt4_7<-m2post$b_Intercept+m2post$b_date0f7+m2post$b_trt0f4+m2post$`b_date0f7:trt0f4` #posterior treatment 4 on date 7


#estimate cumulative emergence. Estimates from poterior are currently in log scale, so take exp to put them on the raw scale.
#Also mulitplied each collection by 7 to remove days, so units change from #/m2/day to #/m2. That's not necessary, but it puts it on
# the same units as benthics

#cumulative emergence after treatments added, units #/m2/d
trt1_tot<-(((exp(trt1_2)*2.78)/4))+(((exp(trt1_3)*2.78)/4))+(((exp(trt1_4)*2.78)/4))+(((exp(trt1_5)*2.78)/4))+(((exp(trt1_6)*2.78)/4))+(((exp(trt1_7)*2.78)/4))
trt2_tot<-(((exp(trt2_2)*2.78)/4))+(((exp(trt2_3)*2.78)/4))+(((exp(trt2_4)*2.78)/4))+(((exp(trt2_5)*2.78)/4))+(((exp(trt2_6)*2.78)/4))+(((exp(trt2_7)*2.78)/4))
trt3_tot<-(((exp(trt3_2)*2.78)/4))+(((exp(trt3_3)*2.78)/4))+(((exp(trt3_4)*2.78)/4))+(((exp(trt3_5)*2.78)/4))+(((exp(trt3_6)*2.78)/4))+(((exp(trt3_7)*2.78)/4))
trt4_tot<-(((exp(trt4_2)*2.78)/4))+(((exp(trt4_3)*2.78)/4))+(((exp(trt4_4)*2.78)/4))+(((exp(trt4_5)*2.78)/4))+(((exp(trt4_6)*2.78)/4))+(((exp(trt4_7)*2.78)/4))

mean(trt1_tot) #mean trt1 
mean(trt2_tot) #mean trt2
mean(trt3_tot) #mean trt3
mean(trt4_tot) #mean trt4

quantile(trt1_tot,probs=c(0.025,0.975)) #quantile trt1
quantile(trt2_tot,probs=c(0.025,0.975)) #quantile trt2
quantile(trt3_tot,probs=c(0.025,0.975)) #quantile trt3
quantile(trt4_tot,probs=c(0.025,0.975)) #quantile trt4

trt= c("control", "atrazine", "selenium", "atrazine + selenium")
avg= c(mean(trt1_tot), mean(trt2_tot), mean(trt3_tot), mean(trt4_tot))
lowq= c(quantile(trt1_tot,probs=c(0.025)), quantile(trt2_tot,probs=c(0.025)), quantile(trt3_tot,probs=c(0.025)), quantile(trt4_tot,probs=c(0.025)))
hiq= c(quantile(trt1_tot,probs=c(0.975)), quantile(trt2_tot,probs=c(0.975)), quantile(trt3_tot,probs=c(0.975)), quantile(trt4_tot,probs=c(0.975)))
cumudf<- data.frame(trt, avg, lowq, hiq)
cumudf
write.csv(cumudf, file="cumulative_emergence071917.csv")

##Cumulative Differences between Treatments##

#Comparing Treatments 1 & 4
diff1_4<-trt1_tot-trt4_tot #difference between treatments estimated at each of the 4000 iterations of the MCMC. Thus, this can be used to describe the probability distribution of differences between these treatements
mean(diff1_4) #There were an average of 90 more insects emerging from trt1 than trt4 per m2 over the 6 week period
quantile(diff1_4,probs=c(0.025,0.975)) #We are 95% certain that the true difference in cumulative emergence between trt1 and 4 was between -28 and 225
sum(diff1_4>0)/4000 #There is a 93% probability that the difference in these treatments was > 0. That's pretty good evidence of a difference.

#Comparing Treatments 1 & 2
diff1_2<-trt1_tot-trt2_tot
mean(diff1_2) #Average 72 more insects emerging from trt1 than trt2 per m2 over 6 week period
quantile(diff1_2, probs=c(0.025,0.975))
sum(diff1_2>0)/4000 #86.6% probability that the difference is >0

#Comparing Treatments 1 & 3
diff1_3<-trt1_tot-trt3_tot
mean(diff1_3) #Average 42 more insects emerging from trt1 than trt3 per m2 over 6 week period
quantile(diff1_3, probs=c(0.025,0.975))
sum(diff1_3>0)/4000 #73% probability that the difference is >0

#Comparing Treatments 2 & 3
diff2_3<-trt3_tot-trt2_tot
mean(diff2_3) #Average 29 more insects emerging from trt3 than trt2 per m2 over the 6 week period
quantile(diff2_3, probs=c(0.025,0.975))
sum(diff2_3>0)/4000 #70% probability that the difference is >0

#Comparing Treatments 2 & 4
diff2_4<-trt2_tot-trt4_tot
mean(diff2_4) #Average 18 more insects emerging from trt2 than trt4 per m2 over 6 week period
quantile(diff2_4, probs=c(0.025,0.975))
sum(diff2_4>0)/4000 #65% probability that the difference is >0

#Comparing Treatments 3 & 4
diff3_4<-trt3_tot-trt4_tot
mean(diff3_4) #Average 48 more insects emerging from trt3 than trt4 per m2 over 6 week period
quantile(diff3_4, probs=c(0.025,0.975))
sum(diff3_4>0)/4000 #80% probability that the difference is >0


##Differences at Week 2##

#Trt 2 & 3
trt2_3_d2<-(((exp(trt3_2)*2.78)/4))-(((exp(trt2_2)*2.78)/4))
mean(trt2_3_d2) # trt2 had an average of 26 fewer insects emerging per m2 per day relative to trt3 on date 2
quantile(trt2_3_d2, probs = c(0.025,0.975))
sum(trt2_3_d2>0)/4000 #We are 93.9% certain that the difference in emergence between trts 2 and 3 on date 2 was greater than zero

##Differences at Week 3##

#Trt 1 & 4 
trt1_4_d3<-(((exp(trt1_3)*2.78)/4))-(((exp(trt4_3)*2.78)/4)) 
mean(trt1_4_d3) #trt4 had an average of 39 fewer insects emerging per m2 per day relative to the control on date 3
quantile(trt1_4_d3,probs=c(0.025,0.975)) 
sum(trt1_4_d3>0)/4000 #We are 95% certain that the difference in emergence between trts 1 and 4 on date 3 was greater than zero.

#Trt 1 & 2
trt1_2_d3<-(((exp(trt1_3)*2.78)/4))-(((exp(trt2_3)*2.78)/4))
mean(trt1_2_d3) # trt2 had an average of 47 fewer insects emerging per m2 per day relative to control on date 3
quantile(trt1_2_d3, probs = c(0.025,0.975))
sum(trt1_2_d3>0)/4000 #We are 98% certain that the difference in emergence between trts 1 and 2 on date 3 was greater than zero

#Trt 1 & 3
trt1_3_d3<-(((exp(trt1_3)*2.78)/4))-(((exp(trt3_3)*2.78)/4))
mean(trt1_3_d3) # trt3 had an average of 41 fewer insects emerging per m2 per day relative to control on date 3
quantile(trt1_3_d3, probs = c(0.025,0.975))
sum(trt1_3_d3>0)/4000 #We are 95% certain that the difference in emergence between trts 1 and 3 on date 3 was greater than zero

#Trt 2 & 3
trt2_3_d3<-(((exp(trt3_3)*2.78)/4))-(((exp(trt2_3)*2.78)/4))
mean(trt2_3_d3) # trt3 had an average of 42 fewer insects emerging per m2 per day relative to control on date 3
quantile(trt2_3_d3, probs = c(0.025,0.975))
sum(trt2_3_d3>0)/4000 #69.5%

#Trt 2 & 4
trt2_4_d3<-(((exp(trt4_3)*2.78)/4))-(((exp(trt2_3)*2.78)/4))
mean(trt2_4_d3) # trt3 had an average of 41 fewer insects emerging per m2 per day relative to control on date 3
quantile(trt2_4_d3, probs = c(0.025,0.975))
sum(trt2_4_d3>0)/4000 #74.3%

#Trt 3 & 4
trt3_4_d3<-(((exp(trt4_3)*2.78)/4))-(((exp(trt3_3)*2.78)/4))
mean(trt3_4_d3) # trt3 had an average of 41 fewer insects emerging per m2 per day relative to control on date 3
quantile(trt3_4_d3, probs = c(0.025,0.975))
sum(trt3_4_d3>0)/4000 #54.5%

##Differences at Week 4##

#Trt 2 & 4
trt2_4_d4<-(((exp(trt2_4)*2.78)/4))-(((exp(trt4_4)*2.78)/4))
mean(trt2_4_d4) #trt4 had an average of 20 fewer insects emerging per m2 per day relative to trt2 on date 4
quantile(trt2_4_d4, probs = c(0.025,0.975)) 
sum(trt2_4_d4>0)/4000 #99.8% certain that the difference in emergence between these treatments on date 4 was greater than zero

#Trt 1 & 4
trt1_4_d4<-(((exp(trt1_4)*2.78)/4))-(((exp(trt4_4)*2.78)/4))
mean(trt1_4_d4) #trt4 had an average of 4 fewer insects emerging per m2 per day relative to control on date 4
quantile(trt1_4_d4, probs = c(0.025,0.975)) #We are 95% certain the true difference in these treatments on date 4 is between -0.35 and 8.6 bugs/m2/day
sum(trt1_4_d4>0)/4000 #97% certain that the difference in emergence between these treatments on date 4 was greater than zero

##Differences at Week 5##

#Trt 1 & 4
trt1_4_d5<-(((exp(trt1_5)*2.78)/4))-(((exp(trt4_5)*2.78)/4))
mean(trt1_4_d5) #trt4 had an average of 6 fewer insects emerging per m2 per day relative to control on date 5
quantile(trt1_4_d5, probs = c(0.025,0.975)) 
sum(trt1_4_d5>0)/4000 #98% certain that the difference in emergence between these treatments on date 5 was greater than zero

#Trt 2 & 4 
trt2_4_d5<-(((exp(trt2_5)*2.78)/4))-(((exp(trt4_5)*2.78)/4))
mean(trt2_4_d5) #trt4 had an average of 7 fewer insects emerging per m2 per day relative to trt2 on date 5
quantile(trt2_4_d5, probs = c(0.025,0.975)) 
sum(trt2_4_d5>0)/4000 #98% certain that the difference in emergence between these treatments on date 5 was greater than zero

#Trt 3 & 4
trt3_4_d5<-(((exp(trt3_5)*2.78)/4))-(((exp(trt4_5)*2.78)/4))
mean(trt3_4_d5) #trt4 had an average of 10 fewer insects emerging per m2 per day relative to trt3 on date 5
quantile(trt3_4_d5, probs = c(0.025,0.975)) 
sum(trt3_4_d5>0)/4000 #99.7% certain that the difference in emergence between these treatments on date 5 was greater than zero

##Differences at Week 6##

#Trt 1 & 2
trt1_2_d6<-(((exp(trt1_6)*2.78)/4))-(((exp(trt2_6)*2.78)/4))
mean(trt1_2_d6) #trt2 had an average of 10 fewer insects emerging per m2 per day relative to control on date 6
quantile(trt1_2_d6, probs = c(0.025,0.975)) 
sum(trt1_2_d6>0)/4000 #99.6% certain that the difference in emergence between these treatments on date 6 was greater than zero

#Trt 1 & 4
trt1_4_d6<-(((exp(trt1_6)*2.78)/4))-(((exp(trt4_6)*2.78)/4))
mean(trt1_4_d6) #trt4 had an average of 10 fewer insects emerging per m2 per day relative to control on date 6
quantile(trt1_4_d6, probs = c(0.025,0.975)) 
sum(trt1_4_d6>0)/4000 #99.5% certain that the difference in emergence between these treatments on date 6 was greater than zero

#Trt 3 & 2
trt3_2_d6<-(((exp(trt3_6)*2.78)/4))-(((exp(trt2_6)*2.78)/4))
mean(trt3_2_d6) #trt2 had an average of 12 fewer insects emerging per m2 per day relative to trt3 on date 6
quantile(trt3_2_d6, probs = c(0.025,0.975)) 
sum(trt3_2_d6>0)/4000 #99.6% certain that the difference in emergence between these treatments on date 6 was greater than zero

#Trt 3 & 4
trt3_4_d6<-(((exp(trt3_6)*2.78)/4))-(((exp(trt4_6)*2.78)/4))
mean(trt3_4_d6) #trt4 had an average of 12 fewer insects emerging per m2 per day relative to trt3 on date 6
quantile(trt3_4_d6, probs = c(0.025,0.975)) 
sum(trt3_4_d6>0)/4000 #99.8% certain that the difference in emergence between these treatments on date 6 was greater than zero

##Differences at Week 7##

#Trt 1 & 4
trt1_4_d7<-(((exp(trt1_7)*2.78)/4))-(((exp(trt4_7)*2.78)/4)) 
mean(trt1_4_d7) #trt4 had an average of 24 fewer insects emerging per m2 per day relative to the control on date 7
quantile(trt1_4_d7,probs=c(0.025,0.975))
sum(trt1_4_d7>0)/4000 #We are 99.1% certain that the difference in emergence between trts 1 and 4 on date 7 was greater than zero.

#Trt 2 & 4
trt2_4_d7<-(((exp(trt2_7)*2.78)/4))-(((exp(trt4_7)*2.78)/4))
mean(trt2_4_d7) #trt4 had an average of 15 fewer insects emerging per m2 per day relative to the trt2 on date 7
quantile(trt2_4_d7,probs=c(0.025,0.975))
sum(trt2_4_d7>0)/4000 #We are 96% certain that the difference in emergence between trts 2 and 3 on date 7 was greater than zero.

#Trt 3 & 4 
trt3_4_d7<-(((exp(trt3_7)*2.78)/4))-(((exp(trt4_7)*2.78)/4)) 
mean(trt3_4_d7) #trt4 had an average of 12 fewer insects emerging per m2 per day relative to the trt3 on date 7
quantile(trt3_4_d7,probs=c(0.025,0.975)) 
sum(trt3_4_d7>0)/4000 #We are 92.7% certain that the difference in emergence between trts 3 and 4 on date 7 was greater than zero.

trt1_3_d7<-(((exp(trt1_7)*2.78)/4))-(((exp(trt3_7)*2.78)/4)) 
mean(trt1_3_d7) #trt4 had an average of 12 fewer insects emerging per m2 per day relative to the trt3 on date 7
quantile(trt1_3_d7,probs=c(0.025,0.975)) 
sum(trt1_3_d7>0)/4000 #We are 84.6% certain that the difference in emergence between trts 1 and 3 on date 7 was greater than zero.


#Differences between selenium date 2 and 3
trt3_d2_d3<- (((exp(trt3_2)*2.78)/4))-(((exp(trt3_3)*2.78)/4)) 
mean(trt3_d2_d3)
quantile(trt3_d2_d3,probs=c(0.025,0.975)) 
sum(trt3_d2_d3>0)/4000 #91.1% prob that diff >0


exp(mean(trt1_3))*2.78/4#108.4
exp(mean(trt2_3))*2.78/4#59.6
exp(mean(trt3_3))*2.78/4#68.2
exp(mean(trt4_3))*2.78/4#70.02
