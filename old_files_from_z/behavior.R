library("lme4")
library("reshape")
library("ggplot2")

filenames <- dir(".","behav.csv")
DATA <- NULL

for (i in 1:length(filenames)) {
  temp <- cbind(SID=substr(filenames[i],1,6),read.csv(filenames[i], header=TRUE))
  temp$RT <- as.numeric(as.character(temp$RT))
  temp <- temp[!(is.na(temp$buttonpress.time)),]
  DATA <- rbind(DATA,temp)
}

# use abs(ineq)
DATA$ineq <- abs(DATA$ineq)

# choose.var = 0 if 10/10 option, 1 otherwise
DATA$choose.var <- 1
DATA[DATA$self.payoff==10 & DATA$other.payoff==10,]$choose.var <- 0

# get variable option payoff
DATA$self.var.payoff <- DATA$self.payoff + DATA$self.foregone - 10
DATA$other.var.payoff <- DATA$other.payoff + DATA$other.foregone - 10
DATA$diff.var.payoff <- abs(DATA$self.var.payoff - DATA$other.var.payoff)

# get variable option ineq
DATA$var.ineq.disadvant <- as.numeric(DATA$other.var.payoff>DATA$self.var.payoff) * (DATA$other.var.payoff-DATA$self.var.payoff)
DATA$var.ineq.advant <- as.numeric(DATA$other.var.payoff<DATA$self.var.payoff) * (DATA$other.var.payoff-DATA$self.var.payoff)

# save
#write.csv(DATA,"DG_behav_all_sub.csv", row.names=FALSE)
#DATA <- read.csv("DG_behav_all_sub.csv", header=TRUE)


# predict choose.var by variable payoff (self + other)
g1.var <- glmer(choose.var ~ self.var.payoff + other.var.payoff + (1|SID), data = DATA, family = "binomial")
summary(g1.var)
#capture.output(summary(g1.var), file = "DG_g1.var_summary.txt")
g1.s02 <- glm(choose.var ~ self.var.payoff + other.var.payoff, data = DATA[DATA$SID=="DG_s02",], family = "binomial")


# predict choose.var by variable self payoff and ineqs
g2.var <- glmer(choose.var ~ self.var.payoff + var.ineq.advant + var.ineq.disadvant + (1|SID) , data = DATA, family = "binomial")
summary(g2.var)
#capture.output(summary(g2.var), file = "DG_g2.var_summary.txt")
g2.s02 <- glm(choose.var ~ self.var.payoff + var.ineq.advant + var.ineq.disadvant, data = DATA[DATA$SID=="DG_s02",], family = "binomial")

# predict choose.var by variable payoff (self + other) in abs(ineq)
g3.var <- glmer(choose.var ~ self.var.payoff + other.var.payoff + abs(self.var.payoff-other.var.payoff) + (1|SID), data = DATA, family = "binomial")
summary(g3.var)
#capture.output(summary(g3.var), file = "DG_g3.var_summary.txt")
g3.s02 <- glm(choose.var ~ self.var.payoff + other.var.payoff + abs(self.var.payoff-other.var.payoff) , data = DATA[DATA$SID=="DG_s02",], family = "binomial")



# attempt plot for g1.var
DATA.subj.self.payoff <- recast(DATA, SID+self.var.payoff~variable, measure.var="choose.var", mean)
DATA.subj.other.payoff <- recast(DATA, SID+other.var.payoff~variable, measure.var="choose.var", mean)
DATA.subj.advant <- recast(DATA, SID+var.ineq.advant~variable, measure.var="choose.var", mean)
DATA.subj.disadvant <- recast(DATA, SID+var.ineq.disadvant~variable, measure.var="choose.var", mean)

names(DATA.subj.self.payoff)[2] <- "amount"
names(DATA.subj.other.payoff)[2] <- "amount"
names(DATA.subj.advant)[2] <- "amount"
names(DATA.subj.disadvant)[2] <- "amount"

DATA.subj.self.payoff$predictor <- "self.payoff"
DATA.subj.other.payoff$predictor <- "other.payoff"
DATA.subj.advant$predictor <- "var.ineq.advant"
DATA.subj.disadvant$predictor <- "var.ineq.disadvant"

DATA.fig <- rbind(DATA.subj.self.payoff,
                  DATA.subj.other.payoff,
                  DATA.subj.advant,
                  DATA.subj.disadvant)

X <- 3*(-25:30)
DATA.fig$amount<- as.numeric(as.character(cut(DATA.fig$amount, breaks=X, labels = X[2:length(X)])))

DATA.fig <- DATA.fig[DATA.fig$predictor %in% c("other.payoff","self.payoff"),]
DATA.s02 <- DATA.fig[DATA.fig$SID=="DG_s02",]

DATA.fig.m <- recast(data.frame(DATA.fig), predictor+amount~variable, measure.var="choose.var", mean)
DATA.fig.m$SE <- recast(data.frame(DATA.fig), predictor+amount~variable, measure.var="choose.var", 
                      function(x)(sd(x)/sqrt(length(x))))[,3]
DATA.fig.m$predictor <- factor(as.character(DATA.fig.m$predictor), levels = c("self.payoff","other.payoff","var.ineq.advant","var.ineq.disadvant"))
g1.fig <- ggplot() + 
  geom_errorbar(data = DATA.fig.m,aes(x=amount, ymin=(choose.var - SE), ymax=(choose.var + SE),color="avg"), 
                position=position_dodge(width=.9), width=0.2) +
  geom_point(data = DATA.fig.m, aes(x=amount, y=choose.var,color="avg")) +
  geom_point(data=DATA.s02, aes(x=amount,y=choose.var,color=SID))+
  facet_grid(.~predictor, scale="free_x") +
  theme_bw(base_size=20) +
  ylab("prop var choice") +
  ylim(c(0,1))+
  xlab("value") +
  theme(strip.background = element_blank())

pdf("DG_prop-var-choice_subjavg.pdf", width = 10, height = 3)
print(g1.fig)
dev.off()
