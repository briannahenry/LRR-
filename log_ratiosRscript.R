#trt2=atrazine; trt3=selenium; trt4=combined#

###log-ratios examples
LR2_3<-log(exp(trt2_3)/exp(trt1_3))
LR3_3<-log(exp(trt3_3)/exp(trt1_3))
summary(LR2_3)
quantile(LR2_3,probs=c(0.025))
mean(LR2_3)

LR2_2<-log(exp(trt2_6)/exp(trt1_6))
LR2_4<-log(exp(trt3_6)/exp(trt1_6))

boxplot(LR2_3,LR3_3,LR2_2,LR2_4,outline=FALSE)
boxplot(LR2_3,LR3_3,LR2_2,outline=FALSE)

#LRR for larvae at time 3 (week 7)
LRL2_3<-log(exp(ltrt2_3)/exp(ltrt1_3))
LRL3_3<-log(exp(ltrt3_3)/exp(ltrt1_3))
LRL4_3<-log(exp(ltrt4_3)/exp(ltrt1_3))

##LRR for cumulative emergence
LRE2<-log(trt2_tot/trt1_tot)
LRE3<-log(trt3_tot/trt1_tot)
LRE4<-log(trt4_tot/trt1_tot)

boxplot(LRL2,LRE2,LRL3,LRE3,LRL4,LRE4,outline=FALSE)

#LRR for emergence at week 3 and week 7
LRE2_3<-log(exp(trt2_3)/exp(trt1_3))
LRE3_3<-log(exp(trt3_3)/exp(trt1_3))
LRE4_3<-log(exp(trt4_3)/exp(trt1_3))

LRE2_7<-log(exp(trt2_7)/exp(trt1_7))
LRE3_7<-log(exp(trt3_7)/exp(trt1_7))
LRE4_7<-log(exp(trt4_7)/exp(trt1_7))

##export to excel
trtLRR= c("larval atrazine_3", "adult atrazine_tot", "adult atrazine_3", "adult atrazine_7",  "larval selenium", "adult selenium_tot", "adult selenium_3", "adult selenium_7", "larval atrazine + selenium", "adult atrazine + selenium_tot", "adult atrazine + selenium_3", "adult atrazine + selenium_7")
LRR= c(mean(LRL2_3), mean(LRE2), mean(LRE2_3), mean(LRE2_7), mean(LRL3_3), mean(LRE3), mean(LRE3_3), mean(LRE3_7), mean(LRL4_3), mean(LRE4), mean(LRE4_3), mean(LRE4_7))
LRRlowq<-c(quantile(LRL2_3,probs=c(0.025)), quantile(LRE2,probs=c(0.025)), quantile(LRE2_3,probs=c(0.025)), quantile(LRE2_7,probs=c(0.025)), quantile(LRL3_3,probs=c(0.025)), quantile(LRE3,probs=c(0.025)), quantile(LRE3_3,probs=c(0.025)), quantile(LRE3_7,probs=c(0.025)), quantile(LRL4_3,probs=c(0.025)), quantile(LRE4,probs=c(0.025)), quantile(LRE4_3,probs=c(0.025)), quantile(LRE4_7,probs=c(0.025)))
LRRhiq<-c(quantile(LRL2_3,probs=c(0.975)), quantile(LRE2,probs=c(0.975)), quantile(LRE2_3,probs=c(0.975)), quantile(LRE2_7,probs=c(0.975)), quantile(LRL3_3,probs=c(0.975)), quantile(LRE3,probs=c(0.975)), quantile(LRE3_3,probs=c(0.975)), quantile(LRE3_7,probs=c(0.975)), quantile(LRL4_3,probs=c(0.975)), quantile(LRE4,probs=c(0.975)), quantile(LRE4_3,probs=c(0.975)), quantile(LRE4_7,probs=c(0.975)))
LRRdf<-data.frame(trtLRR, LRR, LRRlowq, LRRhiq)
write.csv(LRRdf, file="LRR.csv")

#LRR for larvae at time 2 (week 3)
LRL2_2<-log(exp(ltrt2_2)/exp(ltrt1_2))
LRL3_2<-log(exp(ltrt3_2)/exp(ltrt1_2))
LRL4_2<-log(exp(ltrt4_2)/exp(ltrt1_2))

#export larvae and emergence at week 3 to excel
trtLRR= c("larval atrazine_2", "adult atrazine_3", "larval selenium_2","adult selenium_3", "larval atrazine + selenium_2", "adult atrazine + selenium_3")
LRR= c(mean(LRL2_2), mean(LRE2_3), mean(LRL3_2), mean(LRE3_3), mean(LRL4_2),mean(LRE4_3))
LRRlowq<-c(quantile(LRL2_2,probs=c(0.025)),quantile(LRE2_3,probs=c(0.025)), quantile(LRL3_2,probs=c(0.025)), quantile(LRE3_3,probs=c(0.025)),quantile(LRL4_2,probs=c(0.025)), quantile(LRE4_3,probs=c(0.025)))
LRRhiq<-c(quantile(LRL2_2,probs=c(0.975)),quantile(LRE2_3,probs=c(0.975)), quantile(LRL3_2,probs=c(0.975)), quantile(LRE3_3,probs=c(0.975)),quantile(LRL4_2,probs=c(0.975)), quantile(LRE4_3,probs=c(0.975)))
LRR2df<-data.frame(trtLRR, LRR, LRRlowq, LRRhiq)
write.csv(LRR2df, file="LRR2.csv")


###Cumulative emergence additive or not?###
diff2_tot<-trt2_tot-trt1_tot
mean(diff2_tot)
diff3_tot<-trt3_tot-trt1_tot
expect4_tot<-diff2_tot+diff3_tot
mean(expect4_tot)#expected difference from control in total abundance IF trt 2 and 3 are acting independently 

predict4_tot<-trt1_tot+expect4_tot 
mean(predict4_tot)#If trts 2 & 3 were acting independently, then the cumulative bugs for trt 4 should be 649
observed4_tot<-trt4_tot
mean(observed4_tot) #It is actually 742
diff_obs_tot<-observed4_tot-predict4_tot
mean(diff_obs_tot)
sum(diff_obs_tot>0)/4000 #60%
#There is a 60% probability that the observed and predicted are not the same, and it's more likely that 
#the effect is < additive

######observed-expected time 3#####
diff2_3<-exp(trt2_3)-exp(trt1_3)
diff3_3<-exp(trt3_3)-exp(trt1_3)
expect_4<-diff2_3+diff3_3
mean(expect_4)

predict4_3<-exp(trt1_3)+expect_4
observe4_3<-exp(trt4_3)

diff_obs_3<-observe4_3-predict4_3
sum(diff_obs_3>0)/4000 #95% likely trt 2 and 3 not independent 


######observed-expected time 7#####
diff2_7<-exp(trt2_7)-exp(trt1_7)
diff3_7<-exp(trt3_7)-exp(trt1_7)
expect_4_7<-diff2_7+diff3_7
mean(expect_4_7)

predict4_7<-exp(trt1_7)+expect_4_7
observe4_7<-exp(trt4_7)

diff_obs_7<-observe4_7-predict4_7
sum(diff_obs_7>0)/4000 #40% likely that trt 2 and 3 not independent
