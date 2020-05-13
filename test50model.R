sp <- "SEPI_OFF" #MEDITS or SOLEMON code e.g: MERL_MER 

path <- paste("C:/Users/f.masnadi/Desktop/BIOLOGIA della PESCA/CNR/Progetto SEAFAIR/OSR5",sep="")
dir.create(sp) # create carpet with your species name (per me fondamentale avendo più specie)
setwd(paste(path,sp,sep="/")) # set dir to previous carpet 
getwd()

#input <- "your data after Vif" # data (Rdata) after vif saved in previous carpet 
#data.ex <- cbind(readRDS(input),sp)
input <- as.data.frame(input)
class(input)
# fino a qui fate come volete per caricare i vostri dati, era solo un esempio

set.seed(1809) # seleziono un seme
times=50  # numero dei training and test data set
p=0.7 # percentuale per il training 

'------------------Function------------------------------------'
# function to sample train and test == times                  #
f.dyn.tr_ts <- function(data,times){
  classes <- data[,"PA"]
  train_set <- createDataPartition(classes, p = p, list = FALSE,times=times)
  tent <- list()
  mylist.names <- paste("Fold",c(1:times),sep="")
  tent <- vector("list", length(mylist.names))
  names(tent) <- mylist.names
  
  #fun.tr_ts <- function(ro){
  for (i in 1:ncol(train_set)){ 
    tr <- data[train_set[,i],]
    ts <- data[-train_set[,i],]
    
    tr$type <-'Train'
    ts$type <-'Test'
    r <- rbind(tr,ts)
    tent[[i]] <-  r
    #names(tent[[i]]) <- paste("Fold",i,sep="")
    
  }
  return(tent)
}                      #
#
# function to fit DELTA-binomial and to get diagnostic        #
fun.binomial <- function(sb,form){
  gam_train.bin <- gam(formula(form),data=subset(sb,type=='Train'),family=binomial,select=T)
  tes <- data.frame(subset(sb,type=='Test'))
  
  pred.test.bin <- predict(gam_train.bin,newdata=tes,type="response")
  sumry <- summary(gam_train.bin)
  #pred.test_bin=predict(gam_binomial.train,test_ds,type='response')                  # make prediction for test set
  preds.obs=data.frame(pred.test_bin=pred.test.bin,PA=tes$PA)           # data frame of preds vs obs 
  gam1.eval=prediction(preds.obs$pred.test_bin,preds.obs$PA)    # Asses prediction for AUC
  auc.s = attributes(performance(gam1.eval, 'auc'))$y.values[[1]]          # get AUC value
  
  df <- data.frame( AIC=AIC(gam_train.bin),
                    dev.expl=summary(gam_train.bin)$dev.expl,
                    R2.pred = R2(pred.test.bin, tes$PA),
                    #RMSE = RMSE(pred.test.bin,as.data.frame(c[,'PA'])),
                    MAE.pred = MAE(pred.test.bin, tes$PA),
                    AUC.pred=auc.s)
  
  df.list <- list(diag=df,summary=sumry)
  
  my.list <- list(df.list)
  return(my.list)
}                        #
#
# function to fit DELTA-gaussian and to get diagnostic  
# ATTENZIONE!!!!!!!!!!!!!!!!  questa funzione è per i dati
# trasformati con il logaritmo: cambiare prima riga 
# per la trasformazione diversa  (vedi commento interno)      #
fun.gaus_log.eff <- function(sb,form){
  sb2 <- subset(sb,sqrt.kg_km2 >0)# cambiate qui per i dati con radice quadrata se usate quella trasformazione
  c.tr <- data.frame(subset(sb2,type=='Train'))
  c.ts <- data.frame(subset(sb2,type=='Test'))
  gam_train.gaus <- gam(formula(form),data=c.tr)
  pred.test.gaus <- predict(gam_train.gaus,newdata=c.ts,type="response")^2
  
  sumry <- summary(gam_train.gaus)
  
  df <- data.frame( AIC=AIC(gam_train.gaus),
                    BIC=BIC(gam_train.gaus),
                    dev.expl=summary(gam_train.gaus)$dev.expl,
                    R2.pred = R2(pred.test.gaus,c.ts$kg_km2),
                    MAE.pred = MAE(pred.test.gaus, c.ts$kg_km2))
  
  
  df.list <- list(diag=df,summary=sumry,pred= pred.test.gaus)
  
  my.list <- list(df.list)
  return(my.list)
}                    #
#
# function to fit Tweedie and to get diagnostic               #
fun.tw <- function(sb,form){
  gam_train.tw <- gam(formula(form),data=subset(sb,type=='Train'),gamma=1.4,family=tw(),select=T)
  tes <- data.frame(subset(sb,type=='Test'))
  
  pred.test.tw <- predict(gam_train.tw,newdata=tes,type="response")
  sumry <- summary(gam_train.tw)
  #pred.test_bin=predict(gam_binomial.train,test_ds,type='response')                  # make prediction for test set
  preds.obs=data.frame(pred.test_tw=pred.test.tw,kg_km2=tes$kg_km2)           # data frame of preds vs obs 
  #gam1.eval=prediction(preds.obs$pred.test_tw,preds.obs$kg_km2)    # Asses prediction for AUC
  #auc.s = attributes(performance(gam1.eval, 'auc'))$y.values[[1]]          # get AUC value
  
  df <- data.frame( AIC=AIC(gam_train.tw),
                    dev.expl=summary(gam_train.tw)$dev.expl,
                    R2.pred = R2(pred.test.tw, tes$kg_km2),
                    #RMSE = RMSE(pred.test.bin,as.data.frame(c[,'PA'])),
                    MAE.pred = MAE(pred.test.tw, tes$kg_km2))
  #AUC.pred=auc.s)
  
  df.list <- list(diag=df,summary=sumry)
  
  my.list <- list(df.list)
  return(my.list)
}                              #
#
#... function to fit Gaussian and to get diagnostic           #
fun.gaus <- function(sb,form){
  c.tr <- data.frame(subset(sb,type=='Train'))
  c.ts <- data.frame(subset(sb,type=='Test'))
  gam_train.gaus <- gam(formula(form),data=c.tr)
  pred.test.gaus <- predict(gam_train.gaus,newdata=c.ts,type="response")^2
  
  sumry <- summary(gam_train.gaus)
  
  df <- data.frame( AIC=AIC(gam_train.gaus),
                    BIC=BIC(gam_train.gaus),
                    dev.expl=summary(gam_train.gaus)$dev.expl,
                    R2.pred = R2(pred.test.gaus,c.ts$kg_km2),
                    MAE.pred = MAE(pred.test.gaus, c.ts$kg_km2))
  
  
  df.list <- list(diag=df,summary=sumry,pred= pred.test.gaus)
  
  my.list <- list(df.list)
  return(my.list)
}                            #
'-------------------------------------------------------------'

dir.create('2005_2018')
setwd(paste(path,sp,'2005_2018',sep="/")) 
getwd()

# in questa funzione sono presenti delle formule che allego 
# per prima vanno fatte girare queste formule, in base alle vostre variabili selezionate
# con VIF. 
# Probabile per Walter saranno uguali, ma magari Francesco dovrà cambiarle (io per SOLEMON dal vif ho anche POC e ph ma non X)
# lascerei perdere GRAIN visto che è solo GSA 17 

# Formule delta binomial 
formula.eff <- "PA~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(pho)+s(eff)"
formula.1 <- "PA~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(dox.bot)+s(nit)+s(pho)+s(sal)"
formula.2 <- "PA~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(dox.bot)+s(pho)+s(sal)"
formula.3 <- "PA~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(dox.bot)+s(pho)"
formula.4 <- "PA~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(pho)"
formula.5 <- "PA~ s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(pho)"
formula.6 <- "PA~ factor(month)+s(X,Y)+s(year)+s(depth)"

# Formula gaussian (for delta and for all data)
formula.eff.g <- "sqrt.kg_km2~ factor(month)+ s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(nit)+s(eff)"
formula.1.g <- "sqrt.kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(dox.bot)+s(nit)+s(pho)+s(sal)"
formula.2.g <- "sqrt.kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(nit)+s(pho)+s(sal)"
formula.3.g <- "sqrt.kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(nit)+s(pho)"
formula.4.g <- "sqrt.kg_km2~ s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(nit)+s(pho)"
formula.5.g <- "sqrt.kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)"

#  Tweedie 
formula.eff.tw <- "kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(dox.bot)+s(pho)+s(sal)+s(eff)"
formula.1.tw <- "kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(dox.bot)+s(nit)+s(pho)+s(sal)"
formula.2.tw <- "kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(dox.bot)+s(pho)+s(sal)"
formula.3.tw <- "kg_km2~ s(X,Y)+s(year)+s(depth)+s(TMP_bot)+s(dox.bot)+s(pho)+s(sal)"
formula.4.tw <- "kg_km2~ factor(month)+s(X,Y)+s(year)+s(depth)"

# FUNZIONE per 50 trainig/test per il data set dal 2005-2018 (no effort)
# per DELTA-binomial, DELTA-gaussian, Tweedie, Gaussian
mega.fun.noeff <- function(data.ex){
  #.. binomial
  
  trts.all <- f.dyn.tr_ts(data.ex,times)
  var.1 <- sapply(trts.all,fun.binomial,form=formula.1)
  saveRDS(var.1,paste(sp,"deltabin_05_18_var.1","rda",sep="."))
  
  var.2 <- sapply(trts.all,fun.binomial,form=formula.2)
  saveRDS(var.2,paste(sp,"deltabin_05_18_var.2","rda",sep="."))
  
  var.3 <- sapply(trts.all,fun.binomial,form=formula.3)
  saveRDS(var.3,paste(sp,"deltabin_05_18_var.3","rda",sep="."))
  
  var.4 <- sapply(trts.all,fun.binomial,form=formula.4)
  saveRDS(var.4,paste(sp,"deltabin_05_18_var.4","rda",sep="."))
  
  var.5 <- sapply(trts.all,fun.binomial,form=formula.5)
  saveRDS(var.5,paste(sp,"deltabin_05_18_var.5","rda",sep="."))
  
  var.6 <- sapply(trts.all,fun.binomial,form=formula.6)
  saveRDS(var.6,paste(sp,"deltabin_05_18_var.6","rda",sep="."))
  
  df.var.1 <-lapply(var.1, "[[", 1) %>% bind_rows
  df.var.1$mod <- '1'
  
  df.var.2 <- lapply(var.2, "[[", 1) %>% bind_rows
  df.var.2$mod <- '2'
  
  df.var.3 <- lapply(var.3, "[[", 1) %>% bind_rows
  df.var.3$mod <- '3'
  
  df.var.4 <- lapply(var.4, "[[", 1) %>% bind_rows
  df.var.4$mod <- '4'
  
  df.var.5 <- lapply(var.5, "[[", 1) %>% bind_rows 
  df.var.5$mod <- '5'
  
  df.var.6 <- lapply(var.6, "[[", 1) %>% bind_rows 
  df.var.6$mod <- '6'
  
  df.var <- rbind(df.var.1,df.var.2,df.var.3,df.var.4,df.var.5,df.var.6)
  
  melt.df.db <- melt(df.var,id.vars='mod')
  
  melt.df.db2 <- subset(melt.df.db,variable=='AIC'|variable=='dev.expl'|variable=='R2.pred')
  
  p.var.db <- ggplot()+
    geom_boxplot(data=melt.df.db2,aes(x=mod,y=value,fill=mod))+
    facet_wrap(~variable,scale="free_y")+theme(text=element_text(size=18))+
    theme_bw()+ggtitle(paste(sp,"2005:2018 Delta:step binomial train/test"))
  
  #.... delta gaussian ....

  
  var.1.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.1.g)
  saveRDS(var.1.g,paste(sp,"deltaGaus_05_18_var.1.g","rda",sep="."))
  
  var.2.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.2.g)
  saveRDS(var.2.g,paste(sp,"deltaGaus_05_18_var.2.g","rda",sep="."))
  
  var.3.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.3.g)
  saveRDS(var.3.g,paste(sp,"deltaGaus_05_18_var.3.g","rda",sep="."))
  
  var.4.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.4.g)
  saveRDS(var.4.g,paste(sp,"deltaGaus_05_18_var.4.g","rda",sep="."))
  
  var.5.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.5.g)
  saveRDS(var.5.g,paste(sp,"deltaGaus_05_18_var.5.g","rda",sep="."))
  
  df.var.1.g <-lapply(var.1.g, "[[", 1) %>% bind_rows
  df.var.1.g$mod <- '1'
  
  df.var.2.g <- lapply(var.2.g, "[[", 1) %>% bind_rows
  df.var.2.g$mod <- '2'
  
  df.var.3.g <- lapply(var.3.g, "[[", 1) %>% bind_rows
  df.var.3.g$mod <- '3'
  
  df.var.4.g <- lapply(var.4.g, "[[", 1) %>% bind_rows
  df.var.4.g$mod <- '4'
  
  df.var.5.g <- lapply(var.5.g, "[[", 1) %>% bind_rows
  df.var.5.g$mod <- '5'
  
  df.var.dg <- rbind(df.var.1.g,df.var.2.g,df.var.3.g,df.var.4.g,df.var.5.g)
  
  melt.df.dg <- melt(df.var.dg,id.vars='mod')
  
  melt.df.dg2 <- subset(melt.df.dg,variable=='AIC'|variable=="dev.expl"|variable=="R2.pred")
  
  p.var.dg <- ggplot()+
    geom_boxplot(data=melt.df.dg2,aes(x=mod,y=value,fill=mod))+
    facet_wrap(~variable,scale="free_y")+theme(text=element_text(size=18))+
    theme_bw()+ggtitle(paste(sp,"2005:2018 Delta:step Gaussian train/test"))
  
 
  
  ########à end delta gaussian ###############################################
  
  #... Tweedie .... 
  var.1.tw <- sapply(trts.all,fun.tw,form=formula.1.tw)
  saveRDS(var.1.tw,paste(sp,"Tw_05_18_var.1.tw","rda",sep="."))
  
  var.2.tw <- sapply(trts.all,fun.tw,form=formula.2.tw)
  saveRDS(var.2.tw,paste(sp,"Tw_05_18_var.2.tw","rda",sep="."))
  
  var.3.tw <- sapply(trts.all,fun.tw,form=formula.3.tw)
  saveRDS(var.3.tw,paste(sp,"Tw_05_18_var.3.tw","rda",sep="."))
  
  var.4.tw <- sapply(trts.all,fun.tw,form=formula.4.tw)
  saveRDS(var.4.tw,paste(sp,"Tw_05_18_var.4.tw","rda",sep="."))
  
 # df.var.eff.tw <-lapply(var.eff.tw, "[[", 1) %>% bind_rows
#  df.var.eff.tw$mod <- '0'
  
  df.var.1.tw <-lapply(var.1.tw, "[[", 1) %>% bind_rows
  df.var.1.tw$mod <- '1'
  
  df.var.2.tw <- lapply(var.2.tw, "[[", 1) %>% bind_rows
  df.var.2.tw$mod <- '2'
  
  df.var.3.tw <- lapply(var.3.tw, "[[", 1) %>% bind_rows
  df.var.3.tw$mod <- '3'
  
  df.var.4.tw <- lapply(var.4.tw, "[[", 1) %>% bind_rows
  df.var.4.tw$mod <- '4'
  
  df.var.tw <- rbind(df.var.1.tw,df.var.2.tw,df.var.3.tw,df.var.4.tw)
  melt.df.tw <- melt(df.var.tw,id.vars='mod')
  
  melt.df.tw2 <- subset(melt.df.tw,variable=='AIC'|variable=='dev.expl'|variable=='R2.pred')
  
  p.var.tw <- ggplot()+
    geom_boxplot(data=melt.df.tw2,aes(x=mod,y=value,fill=mod))+
    facet_wrap(~variable,scale="free_y")+theme(text=element_text(size=18))+
    theme_bw()+ggtitle(paste(sp,"2005:2018 Result Tweedie train/test"))
  
  ########### end Tweedie ##### 
  
  ############### GAUSSIAN ###########################
  
  var.1.gaus <- sapply(trts.all,fun.gaus,form=formula.1.g)
  saveRDS(var.1.gaus,paste(sp,"Gaus_05_18_var.1.gaus","rda",sep="."))
  
  var.2.gaus <- sapply(trts.all,fun.gaus,form=formula.2.g)
  saveRDS(var.2.gaus,paste(sp,"Gaus_05_18_var.2.gaus","rda",sep="."))
  
  var.3.gaus <- sapply(trts.all,fun.gaus,form=formula.3.g)
  saveRDS(var.3.gaus,paste(sp,"Gaus_05_18_var.3.gaus","rda",sep="."))
  
  var.4.gaus <- sapply(trts.all,fun.gaus,form=formula.4.g)
  saveRDS(var.4.gaus,paste(sp,"Gaus_05_18_var.4.gaus","rda",sep="."))
  
  var.5.gaus <- sapply(trts.all,fun.gaus,form=formula.5.g)
  saveRDS(var.5.gaus,paste(sp,"Gaus_05_18_var.5.gaus","rda",sep="."))
  
  df.var.1.gaus <-lapply(var.1.gaus, "[[", 1) %>% bind_rows
  df.var.1.gaus$mod <- '1'
  
  df.var.2.gaus <- lapply(var.2.gaus, "[[", 1) %>% bind_rows
  df.var.2.gaus$mod <- '2'
  
  df.var.3.gaus <- lapply(var.3.gaus, "[[", 1) %>% bind_rows
  df.var.3.gaus$mod <- '3'
  
  df.var.4.gaus <- lapply(var.4.gaus, "[[", 1) %>% bind_rows
  df.var.4.gaus$mod <- '4'
  
  df.var.5.gaus <- lapply(var.5.gaus, "[[", 1) %>% bind_rows 
  df.var.5.gaus$mod <- '5'
  
  df.var.gaus <- rbind(df.var.1.gaus,df.var.2.gaus,df.var.3.gaus,df.var.4.gaus,df.var.5.gaus)
  
  melt.df.gaus <- melt(df.var.gaus,id.vars='mod')
  
  melt.df.gaus2 <- subset(melt.df.gaus,variable=='AIC'|variable=='dev.expl'|variable=='R2.pred')
  
  
  p.var.gaus <- ggplot()+
    geom_boxplot(data=melt.df.gaus2,aes(x=mod,y=value,fill=mod))+
    facet_wrap(~variable,scale="free_y")+theme(text=element_text(size=18))+
    theme_bw()+ggtitle(paste(sp,"2005:2018 Result Gaussian train/test"))
  
  
  l.fin <- list()
  l.fin[[1]] <- melt.df.db
  l.fin[[2]] <- p.var.db
  l.fin[[3]] <- melt.df.dg 
  l.fin[[4]] <- p.var.dg 
  l.fin[[5]] <- melt.df.tw
  l.fin[[6]] <- p.var.tw
  l.fin[[7]] <- melt.df.gaus
  l.fin[[8]] <- p.var.gaus
  
  names(l.fin) <- c('noEFFmelt.df.db','noEFFboxplot.db','noEFFmelt.df.dg','noEFFboxplot.dg','noEFFmelt.df.tw','noEFFboxplot.tw','noEFFmelt.df.gaus','noEFFboxplot.gaus')
  return(l.fin)
  
}
trts_50.noeff <- mega.fun.noeff(input)

trts_50.noeff$noEFFboxplot.gaus







'<°)))>< <°)))>< <°)))>< <°)))>< <°)))>< <°)))>< <°)))>< <°)))>< <°)))>< <°)))>< '

setwd(paste(path,sp,sep="/")) # set dir to previous carpet 
dir.create('2008_2018')
setwd(paste(path,sp,'2008_2018',sep="/")) 

# FUNZIONE per 50 trainig/test per il data set dal 2008-2018 (con effort)
# per DELTA-binomial, DELTA-gaussian, Tweedie, Gaussian

mega.fun.EFF <- function(data.ex){
  
  data.examin <- (data.ex[data.ex$year>='2008',]%>% dplyr::filter(month > 1)) 
  trts.all <- f.dyn.tr_ts(data.examin,times)
  
  # ... delta binomial ...
  var.eff <- sapply(trts.all,fun.binomial,form=formula.eff)
  saveRDS(var.eff,paste(sp,"deltabin_08_18_var.eff","rda",sep="."))
  
  var.1 <- sapply(trts.all,fun.binomial,form=formula.1)
  saveRDS(var.1,paste(sp,"deltabin_08_18_var.1","rda",sep="."))
  
  var.2 <- sapply(trts.all,fun.binomial,form=formula.2)
  saveRDS(var.2,paste(sp,"deltabin_08_18_var.2","rda",sep="."))
  
  var.3 <- sapply(trts.all,fun.binomial,form=formula.3)
  saveRDS(var.3,paste(sp,"deltabin_08_18_var.3","rda",sep="."))
  
  var.4 <- sapply(trts.all,fun.binomial,form=formula.4)
  saveRDS(var.4,paste(sp,"deltabin_08_18_var.4","rda",sep="."))
  
  var.5 <- sapply(trts.all,fun.binomial,form=formula.5)
  saveRDS(var.5,paste(sp,"deltabin_08_18_var.5","rda",sep="."))
  
  var.6 <- sapply(trts.all,fun.binomial,form=formula.6)
  saveRDS(var.6,paste(sp,"deltabin_08_18_var.6","rda",sep="."))
  
  df.var.eff <-lapply(var.eff, "[[", 1) %>% bind_rows
  df.var.eff$mod <- '0'
  
  df.var.1 <-lapply(var.1, "[[", 1) %>% bind_rows
  df.var.1$mod <- '1'
  
  df.var.2 <- lapply(var.2, "[[", 1) %>% bind_rows
  df.var.2$mod <- '2'
  
  df.var.3 <- lapply(var.3, "[[", 1) %>% bind_rows
  df.var.3$mod <- '3'
  
  df.var.4 <- lapply(var.4, "[[", 1) %>% bind_rows
  df.var.4$mod <- '4'
  
  df.var.5 <- lapply(var.5, "[[", 1) %>% bind_rows 
  df.var.5$mod <- '5'
  
  df.var.6 <- lapply(var.6, "[[", 1) %>% bind_rows 
  df.var.6$mod <- '6'
  
  df.var <- rbind(df.var.eff,df.var.1,df.var.2,df.var.3,df.var.4,df.var.5,df.var.6)
  
  melt.df.db <- melt(df.var,id.vars='mod')
  
  melt.df.db2 <- subset(melt.df.db,variable=='AIC'|variable=='dev.expl'|variable=='R2.pred')
  
  p.var.db <- ggplot()+
    geom_boxplot(data=melt.df.db2,aes(x=mod,y=value,fill=mod))+
    facet_wrap(~variable,scale="free_y")+theme(text=element_text(size=18))+
    theme_bw()+ggtitle(paste(sp,"2008_2018 Delta:step binomial train/test"))
  
  #... delta Gaussian con effort
  
  var.eff.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.eff.g)
  saveRDS(var.eff.g,paste(sp,"08_18_deltaGaus_var.eff.g","rda",sep="."))
  
  var.1.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.1.g)
  saveRDS(var.1.g,paste(sp,"08_18_deltaGaus_var.1.g","rda",sep="."))
  
  var.2.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.2.g)
  saveRDS(var.2.g,paste(sp,"08_18_deltaGaus_var.2.g","rda",sep="."))
  
  var.3.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.3.g)
  saveRDS(var.3.g,paste(sp,"08_18_deltaGaus_var.3.g","rda",sep="."))
  
  var.4.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.4.g)
  saveRDS(var.4.g,paste(sp,"08_18_deltaGaus_var.4.g","rda",sep="."))
  
  var.5.g <- sapply(trts.all,fun.gaus_log.eff,form=formula.5.g)
  saveRDS(var.5.g,paste(sp,"08_18_deltaGaus_var.5.g","rda",sep="."))
  
  df.var.eff.g <-lapply(var.eff.g, "[[", 1) %>% bind_rows
  df.var.eff.g$mod <- '0'
  
  df.var.1.g <-lapply(var.1.g, "[[", 1) %>% bind_rows
  df.var.1.g$mod <- '1'
  
  df.var.2.g <- lapply(var.2.g, "[[", 1) %>% bind_rows
  df.var.2.g$mod <- '2'
  
  df.var.3.g <- lapply(var.3.g, "[[", 1) %>% bind_rows
  df.var.3.g$mod <- '3'
  
  df.var.4.g <- lapply(var.4.g, "[[", 1) %>% bind_rows
  df.var.4.g$mod <- '4'
  
  df.var.5.g <- lapply(var.5.g, "[[", 1) %>% bind_rows
  df.var.5.g$mod <- '5'
  
  df.var.g <- rbind(df.var.eff.g,df.var.1.g,df.var.2.g,df.var.3.g,df.var.4.g,df.var.5.g)
  melt.df.dg <- melt(df.var.g,id.vars='mod')
  
  melt.df.dg2 <- subset(melt.df.dg,variable=='AIC'|variable=="dev.expl"|variable=="R2.pred")
  
  p.var.dg <- ggplot()+
    geom_boxplot(data=melt.df.dg2,aes(x=mod,y=value,fill=mod))+
    facet_wrap(~variable,scale="free_y")+theme(text=element_text(size=18))+
    theme_bw()+ggtitle(paste(sp,"2008_2018 Delta:step Gaussian train/test"))
  ############### END DELTA ####################################################
  
  ###### TWEEDIE   #######################################à
  var.eff.tw <- sapply(trts.all,fun.tw,form=formula.eff.tw)
  saveRDS(var.eff.tw,paste(sp,"TW_08_18_var.eff.tw","rda",sep="."))
  
  var.1.tw <- sapply(trts.all,fun.tw,form=formula.1.tw)
  saveRDS(var.1.tw,paste(sp,"Tw_08_18_var.1.tw","rda",sep="."))
  
  var.2.tw <- sapply(trts.all,fun.tw,form=formula.2.tw)
  saveRDS(var.2.tw,paste(sp,"Tw_08_18_var.2.tw","rda",sep="."))
  
  var.3.tw <- sapply(trts.all,fun.tw,form=formula.3.tw)
  saveRDS(var.3.tw,paste(sp,"Tw_08_18_var.3.tw","rda",sep="."))
  
  var.4.tw <- sapply(trts.all,fun.tw,form=formula.4.tw)
  saveRDS(var.4.tw,paste(sp,"Tw_08_18_var.4.tw","rda",sep="."))
 
  
  df.var.eff.tw <-lapply(var.eff.tw, "[[", 1) %>% bind_rows
  df.var.eff.tw$mod <- '0'
  
  df.var.1.tw <-lapply(var.1.tw, "[[", 1) %>% bind_rows
  df.var.1.tw$mod <- '1'
  
  df.var.2.tw <- lapply(var.2.tw, "[[", 1) %>% bind_rows
  df.var.2.tw$mod <- '2'
  
  df.var.3.tw <- lapply(var.3.tw, "[[", 1) %>% bind_rows
  df.var.3.tw$mod <- '3'
  
  df.var.4.tw <- lapply(var.4.tw, "[[", 1) %>% bind_rows
  df.var.4.tw$mod <- '4'
  
  
  df.var.tw <- rbind(df.var.eff.tw,df.var.1.tw,df.var.2.tw,df.var.3.tw,df.var.4.tw)
  melt.df.tw <- melt(df.var.tw,id.vars='mod')
  
  melt.df.tw2 <- subset(melt.df.tw,variable=='AIC'|variable=='dev.expl'|variable=='R2.pred')
  
  p.var.tw <- ggplot()+
    geom_boxplot(data=melt.df.tw2,aes(x=mod,y=value,fill=mod))+
    facet_wrap(~variable,scale="free_y")+theme(text=element_text(size=18))+
    theme_bw()+ggtitle(paste(sp,"2008:2018 Result Tweedie train/test"))
  
  ############### END TWEEDIE  #################################à
  
  # Gaussian 
  
  var.eff.gaus <- sapply(trts.all,fun.gaus,form=formula.eff.g)
  saveRDS(var.eff.gaus,paste(sp,"Gaus_08_18_var.eff.gaus","rda",sep="."))
  
  var.1.gaus <- sapply(trts.all,fun.gaus,form=formula.1.g)
  saveRDS(var.1.gaus,paste(sp,"Gaus_08_18_var.1.gaus","rda",sep="."))
  
  var.2.gaus <- sapply(trts.all,fun.gaus,form=formula.2.g)
  saveRDS(var.2.gaus,paste(sp,"Gaus_08_18_var.2.gaus","rda",sep="."))
  
  var.3.gaus <- sapply(trts.all,fun.gaus,form=formula.3.g)
  saveRDS(var.3.gaus,paste(sp,"Gaus_08_18_var.3.gaus","rda",sep="."))
  
  var.4.gaus <- sapply(trts.all,fun.gaus,form=formula.4.g)
  saveRDS(var.4.gaus,paste(sp,"Gaus_08_18_var.4.gaus","rda",sep="."))
  
  var.5.gaus <- sapply(trts.all,fun.gaus,form=formula.5.g)
  saveRDS(var.5.gaus,paste(sp,"Gaus_08_18_var.5.gaus","rda",sep="."))
  
  df.var.eff.gaus <-lapply(var.eff.gaus, "[[", 1) %>% bind_rows
  df.var.eff.gaus$mod <- '0'
  
  df.var.1.gaus <-lapply(var.1.gaus, "[[", 1) %>% bind_rows
  df.var.1.gaus$mod <- '1'
  
  df.var.2.gaus <- lapply(var.2.gaus, "[[", 1) %>% bind_rows
  df.var.2.gaus$mod <- '2'
  
  df.var.3.gaus <- lapply(var.3.gaus, "[[", 1) %>% bind_rows
  df.var.3.gaus$mod <- '3'
  
  df.var.4.gaus <- lapply(var.4.gaus, "[[", 1) %>% bind_rows
  df.var.4.gaus$mod <- '4'
  
  df.var.5.gaus <- lapply(var.5.gaus, "[[", 1) %>% bind_rows 
  df.var.5.gaus$mod <- '5'
  
  df.var.gaus <- rbind(df.var.eff.gaus,df.var.1.gaus,df.var.2.gaus,df.var.3.gaus,df.var.4.gaus,df.var.5.gaus)
  
  melt.df.gaus <- melt(df.var.gaus,id.vars='mod')
  
  melt.df.gaus2 <- subset(melt.df.gaus,variable=='AIC'|variable=='dev.expl'|variable=='R2.pred')
  
  
  p.var.gaus <- ggplot()+
    geom_boxplot(data=melt.df.gaus2,aes(x=mod,y=value,fill=mod))+
    facet_wrap(~variable,scale="free_y")+theme(text=element_text(size=18))+
    theme_bw()+ggtitle(paste(sp,"2008:2018 Result Gaussian train/test"))
  
  
  
  l.fin <- list()
  l.fin[[1]] <- melt.df.db
  l.fin[[2]] <- p.var.db
  l.fin[[3]] <- melt.df.dg
  l.fin[[4]] <- p.var.dg
  l.fin[[5]] <- melt.df.tw
  l.fin[[6]] <- p.var.tw 
  l.fin[[7]] <- melt.df.gaus
  l.fin[[8]] <- p.var.gaus
  names(l.fin) <- c('EFFmelt.df.db','EFFboxplot.db','EFFmelt.df.dg','EFFboxplot.dg','EFFmelt.df.tw','Effboxplot.tw','Effmelt.df.gaus','Effboxplot.gaus')
  return(l.fin)
  
} # end function
trts_50.eff <- mega.fun.EFF(input)

trts_50.eff$EFFboxplot.dg
# END
'<°)))>< <°)))>< <°)))>< <°)))>< <°)))>< <°)))>< <°)))>< <°)))>< <°)))>< <°)))>< '



