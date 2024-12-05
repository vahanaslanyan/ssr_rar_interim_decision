set.seed(919)

source("bayesian_sse_rar.R")
df1_aux<-df1
df2_aux<-df2
df3_aux<-df3
df4_aux<-df4
#data cleaning
df1<-df1%>%arrange(active_treatment)
df1$active_treatment<-as.factor(df1$active_treatment)
df2<-df1
df1$sd_treatment1<-NULL
df2$sd_treatment2<-NULL
df1$treatment1<-NULL
df1$treatment2<-NULL

df1 %>%select(active_treatment,futility,efficacy,
              stopped_for_efficacy_1,
              stopped_for_efficacy_2)%>% pivot_longer(cols=c('stopped_for_efficacy_1','stopped_for_efficacy_2'),
                    names_to='stopped',
                    values_to='points')->df2

df1 %>%select(active_treatment,futility,efficacy,
              stopped_for_futility_1,
              stopped_for_futility_2)%>% 
  pivot_longer(cols=c('stopped_for_futility_1','stopped_for_futility_2'),
                                                      names_to='stopped',
                                                      values_to='points')->df3

df4<-df2%>%bind_rows(df3)
df4$points<-df4$points*100
df4$Run<-rep(c("25%","50%"),56)
df4$Run<-as.factor(df4$Run)
df4$Design<-paste0(df4$efficacy,"-",df4$futility)
df4$Reason<-gsub("stopped_for_","",df4$stopped)
df4$Reason<-gsub("_1","",df4$Reason)
df4$Reason<-gsub("_2","",df4$Reason)

#panel A
p1<-ggplot(df4,aes(x=Run,y=points,color=Design))+
  facet_wrap(~active_treatment+Reason,strip.position = "bottom",ncol=14)+
  geom_point(position = "jitter",size=4)+
  theme_minimal()+
  ylab("% terminated")+xlab("Interim Analysis")+
  theme(legend.position = "bottom")+theme(text = element_text(size=30))+
 scale_color_manual(values=c("#8d230f", "#1e434c","#f0810f","#a1d6e2"))+
  theme(panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))
p1



#repeated from above, same as for all other conditions
df2<-df2_aux
df2<-df2%>%arrange(active_treatment)
df2$active_treatment<-as.factor(df2$active_treatment)

df2$sd_treatment1<-NULL
df2$sd_treatment2<-NULL
df2$treatment1<-NULL
df2$treatment2<-NULL

df2 %>%select(active_treatment,futility,efficacy,
              stopped_for_efficacy_1,
              stopped_for_efficacy_2)%>% pivot_longer(cols=c('stopped_for_efficacy_1','stopped_for_efficacy_2'),
                                                      names_to='stopped',
                                                      values_to='points')->df3

df2 %>%select(active_treatment,futility,efficacy,
              stopped_for_futility_1,
              stopped_for_futility_2)%>% 
  pivot_longer(cols=c('stopped_for_futility_1','stopped_for_futility_2'),
               names_to='stopped',
               values_to='points')->df4

df5<-df3%>%bind_rows(df4)
df5$points<-df5$points*100
df5$Run<-rep(c("25%","50%"),56)
df5$Run<-as.factor(df5$Run)
df5$Design<-paste0(df5$efficacy,"-",df5$futility)
df5$Reason<-gsub("stopped_for_","",df5$stopped)
df5$Reason<-gsub("_1","",df5$Reason)
df5$Reason<-gsub("_2","",df5$Reason)

#panel B
p2<-ggplot(df5,aes(x=Run,y=points,color=Design))+
  facet_wrap(~active_treatment+Reason,strip.position = "bottom",ncol=14)+
  geom_point(position = "jitter",size=4)+
  theme_minimal()+
  ylab("% terminated")+xlab("Interim Analysis")+
  theme(legend.position = "bottom")+theme(text = element_text(size=30))+
  scale_color_manual(values=c("#8d230f", "#1e434c","#f0810f","#a1d6e2"))+
  theme(panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))
p2







df2<-df3_aux
df2<-df2%>%arrange(active_treatment)
df2$active_treatment<-as.factor(df2$active_treatment)

df2$sd_treatment1<-NULL
df2$sd_treatment2<-NULL
df2$treatment1<-NULL
df2$treatment2<-NULL

df2 %>%select(active_treatment,futility,efficacy,
              stopped_for_efficacy_1,
              stopped_for_efficacy_2)%>% pivot_longer(cols=c('stopped_for_efficacy_1','stopped_for_efficacy_2'),
                                                      names_to='stopped',
                                                      values_to='points')->df3

df2 %>%select(active_treatment,futility,efficacy,
              stopped_for_futility_1,
              stopped_for_futility_2)%>% 
  pivot_longer(cols=c('stopped_for_futility_1','stopped_for_futility_2'),
               names_to='stopped',
               values_to='points')->df4

df5<-df3%>%bind_rows(df4)
df5$points<-df5$points*100
df5$Run<-rep(c("25%","50%"),56)
df5$Run<-as.factor(df5$Run)
df5$Design<-paste0(df5$efficacy,"-",df5$futility)
df5$Reason<-gsub("stopped_for_","",df5$stopped)
df5$Reason<-gsub("_1","",df5$Reason)
df5$Reason<-gsub("_2","",df5$Reason)

#panel C
p3<-ggplot(df5,aes(x=Run,y=points,color=Design))+
  facet_wrap(~active_treatment+Reason,strip.position = "bottom",ncol=14)+
  geom_point(position = "jitter",size=4)+
  theme_minimal()+
  ylab("% terminated")+xlab("Interim Analysis")+
  theme(legend.position = "bottom")+theme(text = element_text(size=30))+
  scale_color_manual(values=c("#8d230f", "#1e434c","#f0810f","#a1d6e2"))+
  theme(panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))
p3




df2<-df4_aux
df2<-df2%>%arrange(active_treatment)
df2$active_treatment<-as.factor(df2$active_treatment)

df2$sd_treatment1<-NULL
df2$sd_treatment2<-NULL
df2$treatment1<-NULL
df2$treatment2<-NULL

df2 %>%select(active_treatment,futility,efficacy,
              stopped_for_efficacy_1,
              stopped_for_efficacy_2)%>% pivot_longer(cols=c('stopped_for_efficacy_1','stopped_for_efficacy_2'),
                                                      names_to='stopped',
                                                      values_to='points')->df3

df2 %>%select(active_treatment,futility,efficacy,
              stopped_for_futility_1,
              stopped_for_futility_2)%>% 
  pivot_longer(cols=c('stopped_for_futility_1','stopped_for_futility_2'),
               names_to='stopped',
               values_to='points')->df4

df5<-df3%>%bind_rows(df4)
df5$points<-df5$points*100
df5$Run<-rep(c("25%","50%"),56)
df5$Run<-as.factor(df5$Run)
df5$Design<-paste0(df5$efficacy,"-",df5$futility)
df5$Reason<-gsub("stopped_for_","",df5$stopped)
df5$Reason<-gsub("_1","",df5$Reason)
df5$Reason<-gsub("_2","",df5$Reason)

#Panel D
p4<-ggplot(df5,aes(x=Run,y=points,color=Design))+
  facet_wrap(~active_treatment+Reason,strip.position = "bottom",ncol=14)+
  geom_point(position = "jitter",size=4)+
  theme_minimal()+
  ylab("% terminated")+xlab("Interim Analysis")+
  theme(legend.position = "bottom")+theme(text = element_text(size=30))+
  scale_color_manual(values=c("#8d230f", "#1e434c","#f0810f","#a1d6e2"))+
  theme(panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))
p4







library(ggpubr)
#Figure 1
p5<-ggarrange(p1,p2,p3,p4, nrow=4,common.legend = T, legend="bottom",labels = "AUTO",
          font.label=list(color="black",size=36))

p5


#ggsave("percent_stop.png",p5,dpi=600,width = 22,height = 28)



df2<-df1_aux

df2<-df2%>%arrange(active_treatment)
df2$active_treatment<-as.factor(df2$active_treatment)
df2$stopped_for_efficacy_1<-NULL
df2$stopped_for_efficacy_2<-NULL
df2$stopped_for_futility_1<-NULL
df2$stopped_for_futility_2<-NULL
df2$Design<-paste0(df2$efficacy,"-",df2$futility)
df2 %>%select(active_treatment,treatment1,
              treatment2,Design)->df3
df3%>%pivot_longer(cols=c('treatment1','treatment2'),
                   names_to='treatment',
                   values_to='sample_size')->df4


p6<-ggplot(df4,aes(x=treatment,y=sample_size,fill=Design))+
  facet_wrap(~active_treatment,strip.position = "bottom",ncol=14)+
  geom_col(position="dodge")+
  theme_minimal()+
  ylab("Sample Size")+xlab("Treatment")+
  scale_fill_manual(values=c("#8d230f", "#1e434c","#f0810f","#a1d6e2"))+
  theme(legend.position = "bottom")+theme(text = element_text(size=30))+
  theme(panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))+
  scale_x_discrete(labels= c("Placebo", "Active"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))



df2<-df2_aux

df2<-df2%>%arrange(active_treatment)
df2$active_treatment<-as.factor(df2$active_treatment)
df2$stopped_for_efficacy_1<-NULL
df2$stopped_for_efficacy_2<-NULL
df2$stopped_for_futility_1<-NULL
df2$stopped_for_futility_2<-NULL
df2$Design<-paste0(df2$efficacy,"-",df2$futility)
df2 %>%select(active_treatment,treatment1,
              treatment2,Design)->df3
df3%>%pivot_longer(cols=c('treatment1','treatment2'),
                   names_to='treatment',
                   values_to='sample_size')->df4


p7<-ggplot(df4,aes(x=treatment,y=sample_size,fill=Design))+
  facet_wrap(~active_treatment,strip.position = "bottom",ncol=14)+
  geom_col(position="dodge")+
  theme_minimal()+
  ylab("Sample Size")+xlab("Treatment")+
  scale_fill_manual(values=c("#8d230f", "#1e434c","#f0810f","#a1d6e2"))+
  theme(legend.position = "bottom")+theme(text = element_text(size=30))+
  theme(panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))+
  scale_x_discrete(labels= c("Placebo", "Active"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


df2<-df3_aux

df2<-df2%>%arrange(active_treatment)
df2$active_treatment<-as.factor(df2$active_treatment)
df2$stopped_for_efficacy_1<-NULL
df2$stopped_for_efficacy_2<-NULL
df2$stopped_for_futility_1<-NULL
df2$stopped_for_futility_2<-NULL
df2$Design<-paste0(df2$efficacy,"-",df2$futility)
df2 %>%select(active_treatment,treatment1,
              treatment2,Design)->df3
df3%>%pivot_longer(cols=c('treatment1','treatment2'),
                   names_to='treatment',
                   values_to='sample_size')->df4


p8<-ggplot(df4,aes(x=treatment,y=sample_size,fill=Design))+
  facet_wrap(~active_treatment,strip.position = "bottom",ncol=14)+
  geom_col(position="dodge")+
  theme_minimal()+
  ylab("Sample Size")+xlab("Treatment")+
  scale_fill_manual(values=c("#8d230f", "#1e434c","#f0810f","#a1d6e2"))+
  theme(legend.position = "bottom")+theme(text = element_text(size=30))+
  theme(panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))+
  scale_x_discrete(labels= c("Placebo", "Active"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


df2<-df4_aux

df2<-df2%>%arrange(active_treatment)
df2$active_treatment<-as.factor(df2$active_treatment)
df2$stopped_for_efficacy_1<-NULL
df2$stopped_for_efficacy_2<-NULL
df2$stopped_for_futility_1<-NULL
df2$stopped_for_futility_2<-NULL
df2$Design<-paste0(df2$efficacy,"-",df2$futility)
df2 %>%select(active_treatment,treatment1,
              treatment2,Design)->df3
df3%>%pivot_longer(cols=c('treatment1','treatment2'),
                   names_to='treatment',
                   values_to='sample_size')->df4


p9<-ggplot(df4,aes(x=treatment,y=sample_size,fill=Design))+
  facet_wrap(~active_treatment,strip.position = "bottom",ncol=14)+
  geom_col(position="dodge")+
  theme_minimal()+
  ylab("Sample Size")+xlab("Treatment")+
  scale_fill_manual(values=c("#8d230f", "#1e434c","#f0810f","#a1d6e2"))+
  theme(legend.position = "bottom")+theme(text = element_text(size=30))+
  theme(panel.spacing = unit(.05, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1))+
  theme(plot.background = element_rect(fill = 'white', colour = 'white'))+
  scale_x_discrete(labels= c("Placebo", "Active"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))



library(ggpubr)
#Figure 2
p10<-ggarrange(p6,p7,p8,p9, nrow=4,common.legend = T, legend="bottom",labels = "AUTO",
              font.label=list(color="black",size=36))

p10


#ggsave("sample sizes.png",p10,dpi=600,width = 23,height = 29)
