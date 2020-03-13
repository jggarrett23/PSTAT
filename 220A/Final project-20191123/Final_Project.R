library(faraway)
library(tidyverse)
library(ggplot2)
library(psych)
library(GGally)
library(MASS)
library(car)
library(lava)
library(boot)
library(latex2exp)
library(plyr)
library(cowplot)
library(glmnet)
plot_dir <- '/Users/owner/Desktop/PSTAT_220A/Final project-20191123/Figures'

#====================================== Description of Data Variables ========================================
property.Data <- read.csv("/Users/owner/Desktop/PSTAT_220A/Final project-20191123/Final_Data.csv",
                          sep=",", header=T)


#1. size: size of the property (in square meters).
#2. age: age of the property (in years).
#3. dc: distance (in km) from the property to the city center.
#4. dt: distance (in km) from the property to a toxic waste disposal site.
#5. price: the listed price of the property, in thousands of dollars.


#Investigate how listed price depends on other variables.

descript.stats <- summary(property.Data)

#============================================ Plot Relationships =============================================

cov.boxs <- ggplot(stack(property.Data), aes(x = ind, y = values, fill = ind)) +
  geom_boxplot(outlier.color = 'red') + scale_x_discrete(limits = c("size", "age", "dc", "dt")) +
  theme(legend.position = "none", axis.title = element_blank(), text = element_text(size = 15))

resp.boxs <- ggplot(stack(property.Data), aes(x = ind, y = values)) +
  geom_boxplot(outlier.color = 'red') + scale_x_discrete(limits = c("price")) + theme(axis.title=element_blank(),
                                                                                      text=element_text(size = 15))


box.grid <- plot_grid(cov.boxs,resp.boxs)
ggsave('BoxPlots.jpg',plot=box.grid,path=plot_dir, width = 7, height = 6, units = "cm")


panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
 
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(property.Data, diag.panel = panel.hist, lower.panel=panel.smooth,
      upper.panel = panel.cor, cex.labels=1,font.labels = 2)

my_fn <- function(data, mapping, pts=list(), smt=list(), ...){
  ggplot(data = data, mapping = mapping, ...) + 
    do.call(geom_point, pts) +
    do.call(geom_smooth, smt) 
}

pair.Plot <- ggpairs(property.Data, progress=F,
        lower = list(continuous = 
                       wrap(my_fn,
                            pts=list(size=1, colour="black"), 
                            smt=list(method="loess", se=F, size=1, colour="red"))),
        upper = list(continuous = wrap("cor", color = "blue", size = 5)))
ggsave('Explore_Pairs.jpg',plot=pair.Plot,path=plot_dir, width = 16, height= 11, units = 'cm')

#========================================= Data Cleaning =========================================================

outliers <- boxplot(property.Data, plot=F)$out

#clean.propData <- property.Data[-which(property.Data$price %in% outliers),]
clean.propData <- property.Data[-which(property.Data$age %in% outliers),]

ggpairs(clean.propData, progress=F,
        lower = list(continuous = 
                       wrap(my_fn,
                            pts=list(size=1, colour="black"), 
                            smt=list(method="loess", se=F, size=1, colour="red"))),
        upper = list(continuous = wrap("cor", color = "blue", size = 5)))



#========================================= Fit Models =========================================================

summary(full_fit <- lm(price ~., data = property.Data))

##partial fits

#create T/F vector for presence or absence of regressors
cov <- names(property.Data[c(1:4)])
regMat <- expand.grid(c(T,F), c(T,F), c(T,F), c(T,F))
regMat <- regMat[c(2:15),] #remove first and last row
names(regMat) <- cov

models.List <- apply(regMat, 1, 
                     function(x) as.formula(paste("price ~",paste(c(cov[x]), collapse = '+'))))

models.results <- lapply(models.List,function(x) lm(x, data=property.Data))

dfCoefNum   <- ldply(models.results, function(x) as.data.frame(
  t(coef(x))))
dfStdErrors <- ldply(models.results, function(x) as.data.frame(
  t(coef(summary(x))[, "Std. Error"])))
dftValues   <- ldply(models.results, function(x) as.data.frame(
  t(coef(summary(x))[, "t value"])))
dfpValues   <- ldply(models.results, function(x) as.data.frame(
  t(coef(summary(x))[, "Pr(>|t|)"]))) 

# rename DFs so we know what the column contains
names(dfStdErrors) <- paste("se", names(dfStdErrors), sep=".")
names(dftValues) <- paste("t", names(dftValues), sep=".")
names(dfpValues) <- paste("p", names(dfpValues), sep=".")

# p-value for overall model fit
calcPval <- function(x){
  fstat <- summary(x)$fstatistic
  pVal <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
  return(pVal)
}

# Before creating ONE data frame with all important entries,
# we need to compute some more indices 
NoOfCoef <- unlist(apply(regMat, 1, sum))
R2       <- unlist(lapply(models.results, function(x)
  summary(x)$r.squared))
adjR2    <- unlist(lapply(models.results, function(x)
  summary(x)$adj.r.squared))
RMSE     <- unlist(lapply(models.results, function(x)
  summary(x)$sigma))
fstats   <- unlist(lapply(models.results, calcPval))

results <- data.frame( model = as.character(models.List),
                       NoOfCoef = NoOfCoef,
                       dfCoefNum,
                       dfStdErrors,
                       dftValues,
                       dfpValues,
                       R2 = R2,
                       adjR2 = adjR2,
                       RMSE = RMSE,
                       pF = fstats  )

top5.RSME <- head(results[order(results$RMSE),])

# ========================================= Best Model Diagnostics ===========================================
#extracts the model with minimum RSME
summary(best_fit.1 <- models.results[[rownames(top5.RSME)[1]]])
r <- residuals(trans.best_fit.1)
yh <- trans.best_fit.1$fitted.values
s_r <- rstandard(trans.best_fit.1)

jpeg(paste(plot_dir,"residuals_plots.jpg",sep="/"), width = 352, height = 323, units='px')
par(mfrow=c(2,2))
#jpeg(paste(plot_dir,"residuals_plot1.jpg",sep="/"), width = 352, height = 323, units='px')
plot(yh,s_r, ylab = "Standardized Residuals", xlab = "Fitted Values", 
     main = "Standardized Residuals\n vs Fitted Values")
lines(lowess(yh,s_r),lty=2,col="red")
abline(h=0,lty=3,col='gray')
#dev.off()

#jpeg(paste(plot_dir,"residuals_plot2.jpg",sep="/"), width = 352, height = 323, units='px')
qqnorm(s_r,main="Normal Q-Q Plot", ylab="Standardized Residuals")
qqline(s_r)
#dev.off()

#jpeg(paste(plot_dir,"residuals_plot3.jpg",sep="/"), width = 352, height = 323, units='px')
plot(r[-length(r)],r[-1],
     xlab=expression(hat(epsilon)[i]),
     ylab=expression(hat(epsilon)[i+1]),
     main = "Independence of Residuals", cex.lab = 1.4)
lines(lowess(r[-length(r)],r[-1]), col="red", lty=2)
#dev.off()

#jpeg(paste(plot_dir,"residuals_plot4.jpg",sep="/"), width = 352, height = 323, units='px')
influencePlot(best_fit.1, main = "Influence Plot of Residuals")
dev.off()

# formula: price ~ size + age + dt
# Normality seems fine. Constant variance is a bit worrisome. Outliers/Leverage fine
ncvTest(best_fit.1)
par(mfrow=c(1,1))
influencePlot(best_fit.1)
outlierTest(best_fit.1)

#partial residuals
plot_partials <- function (current_var,model,df) {
  pr <- residuals(model)+coef(model)[[current_var]]*df[[current_var]]
  ggplot(df, aes(x=df[[current_var]], y=pr)) + 
    geom_point(shape = 1) + labs(x=current_var,y='Partial Residuals',title=paste(current_var,'Partial Residuals', sep=' ')) + 
    geom_smooth(method='lm', se = F) +
    geom_smooth(method='loess',se=F, color ='red', linetype='dashed') +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
}


partial.plots <- lapply(names(best_fit.1$coefficients)[c(2:4)],
                        plot_partials,model=best_fit.1,df=property.Data)

grid.partials <- plot_grid(partial.plots[[1]],partial.plots[[2]],partial.plots[[3]])

ggsave('Best_Fit_1_PR.jpg',plot=grid.partials,path=plot_dir)
# age and dt have non linear partial residual relationship, suggesting that transformation might be necessary

# added variables
d.size <- residuals(lm(price ~ dt + age, data=property.Data))
m.size <- residuals(lm(size ~ dt + age, data=property.Data))

size_adv <- ggplot(property.Data, aes(x=m.size, y=d.size)) + 
  geom_point(shape = 1) + labs(x="Size residuals",y='Price residuals',title="Added variable plot for Size") + 
  geom_smooth(method='lm', se = F) +
  geom_smooth(method='loess',se=F, color ='red', linetype='dashed') +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

d.dt <- residuals(lm(price ~ size + age, data=property.Data))
m.dt <- residuals(lm(dt ~ size + age, data=property.Data))

dt_adv <- ggplot(property.Data, aes(x=m.dt, y=d.dt)) + 
  geom_point(shape = 1) + labs(x="Dt residuals",y='Price residuals',title="Added variable plot for Dt") + 
  geom_smooth(method='lm', se = F) +
  geom_smooth(method='loess',se=F, color ='red', linetype='dashed') +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

d.age<- residuals(lm(price ~ size + dt, data=property.Data))
m.age <- residuals(lm(age ~ size + dt, data=property.Data))

age_adv <- ggplot(property.Data, aes(x=m.age, y=d.age)) + 
  geom_point(shape = 1) + labs(x="Age residuals",y='Price residuals',title="Added variable plot for Age") + 
  geom_smooth(method='lm', se = F) +
  geom_smooth(method='loess',se=F, color ='red', linetype='dashed') +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))


grid.addv <- plot_grid(size_adv,dt_adv,age_adv)

# =============================== Control for Collinearity Between Variables ================================
vif(best_fit.1)
vif(full_fit)

# variance inflation not an issue in both models

# ========================================== Step-wise Regression ===========================================

step(lm(price~1, data=property.Data), scope=list(upper=formula(full_fit)), direction="forward")

# =============================================== Transformations ===========================================
boxcox(full_fit,plotit = T, lambda = seq(-0.4,1.2,.2)) #zero in confidence interval, do log transformation
boxcox(lm(price~size+dt+age,data=property.Data),plotit=T) #zero in confidence interval, do log transformation

summary(trans.full_fit <- update(full_fit,I(log(price)) ~ .))
summary(trans.best_fit.1 <- update(best_fit.1, I(log(price)) ~ .))

step(trans.full_fit, direction = 'both')
step(trans.best_fit.1, direction = 'both')


# =================================== Ridge Regression & Lasso ===============================================

#center and scale to mean of 0 and sd of 1
scaled.df <- as.data.frame(scale(property.Data))
summary(full_fit.scaled <- lm(price ~ . , data = scaled.df))
summary(best_fit.1_scaled <- update(full_fit.scaled, . ~ .-dc))

X <- model.matrix(full_fit.scaled)[,-1]
fit.ridge <- glmnet(X, scaled.df$price, lambda.min=0, nlambda=101, alpha=0)
plot(fit.ridge, xvar="lambda", xlim=c(-5,7)) 
text(-1,coef(fit.ridge)[-1,length(fit.ridge$lambda)],labels=colnames(X),cex=0.6) 
fit.ridge.cv <- cv.glmnet(X, scaled.df$price, lambda.min=0, nlambda=101, alpha=0) 
abline(v=log(fit.ridge.cv$lambda.min), col="red")
mtext("CV estimate", side=1, at=log(fit.ridge.cv$lambda.min), cex=.6)
plot(fit.ridge.cv)

set.seed(294)
fit.lasso <- glmnet(X, scaled.df$price, lambda.min=0, nlambda=101, alpha=1)
jpeg(paste(plot_dir,"Lasso_plot.jpg",sep="/"), width=584, height = 550, units = "px")
plot(fit.lasso, xvar="lambda", xlim=c(-6,1), ylim = c(-0.55,0.62)) 
text(-5.5,coef(fit.lasso)[-1,length(fit.lasso$lambda)]-0.05,labels=colnames(X),cex=1.1) 
mtext("CV estimate", side=1, at=log(fit.lasso.cv$lambda.min), cex=.8)
abline(v=log(fit.lasso.cv$lambda.min), col="red")
dev.off()

fit.lasso.cv <- cv.glmnet(X, scaled.df$price, lambda.min=0, nlambda=101) 
plot(fit.lasso.cv)

# ================================== Normal vs Transformed Model  ===========================================
final_models.list <- list(full_fit,trans.full_fit,best_fit.1,trans.best_fit.1)
model.AIC <- unlist(lapply(final_models.list,AIC))
model.BIC <- unlist(lapply(final_models.list,BIC))

model.performance.sum <- data.frame("Formula"=unlist(lapply(final_models.list,function(x)
                                                            Reduce(paste,deparse(formula(x))))),
                                    "R2"=unlist(lapply(final_models.list,function(x)
                                                        summary(x)$r.squared)), 
                                    "AdjR"=unlist(lapply(final_models.list,function(x)
                                      summary(x)$adj.r.squared)),
                                    "RSME" = unlist(lapply(final_models.list, function(x)
                                      summary(x)$sigma)),
                                    "pF"=unlist(lapply(final_models.list, calcPval)),
                                    "AIC" = model.AIC,"BIC"=model.BIC)

# ========================================= Final Diagnostics =========================================
#Note that models with interaction terms for the full fit may have lower AIC, but are more complex & 
#harder to interpret

par(mfrow=c(2,2))
plot(trans.best_fit.1, which = c(1:3,5))

outlierTest(trans.best_fit.1)

partials.2 <- lapply(names(trans.best_fit.1$coefficients)[c(2:4)],
       plot_partials,model=trans.best_fit.1,df=property.Data)

grid.partials2  <- plot_grid(partials.2[[1]],partials.2[[2]],partials.2[[3]])


# ========================================= Outliers Removed =========================================
step(lm(price ~ ., data = clean.propData), direction = "both")

clean.best_fit <- lm(price ~ size + dt, data = clean.propData)

boxcox(clean.best_fit)

trans.cleanBest <- lm(I(log(price)) ~ size + dt, data = clean.propData)

partials.3 <- lapply(names(trans.cleanBest$coefficients)[c(2:3)],
                     plot_partials,model=trans.cleanBest,df=clean.propData)

grid.partials2  <- plot_grid(partials.3[[1]],partials.3[[2]])


# =========================================== Predictions ===========================================
confidenceEllipse(best_fit.1, level=.95,
                  main = "95% Confidence Region for\n Publication and Experience Estimates")
abline(v=confint(best_fit.1, level=.95)[2,], lty=2)
abline(h=confint(best_fit.1, level=.95)[3,], lty=2)
points(0,0, col='red', cex=3, pch=13)
