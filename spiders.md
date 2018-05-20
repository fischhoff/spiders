spiders, Met52, and ticks
================
Ilya
5/19/2018

#### Study design

##### Study question: Are Met52 (fungal biopesticide) and brush-legged wolf spiders (Schizocosa ocreata) compatible biocontrol agents against ticks?

##### Decision: How many microcosms are necessary to detect interaction between Met52 and spider?

##### Approach: bootstrap power analysis using data from previous experiment and assumed effects of Met52 and of Met52*spider interaction. Power analysis: Use data from previous experiment (conducted by collaborator J. Burtis), on survival of ticks in microcosms with and without S. ocreata, to determine sample size needed for current experiment. For range of sample sizes, simulate four treatments: 1) spider + biopesticide control (H2O), 2) spider + Met52, 3) no spider + H2O, 4) no spider + Met52. For biopesticide control microcosms, sample, with replacement, from dataset of tick survival in microcosms with and without spider. For microcosms with biopesticide and no spider, simulate 50% reduction in survival due to Met52 after sampling from data on survival in no-spider microcosms. For microcosms with biopesticide and spider, assign spider effect based on random sample of microcosms with and without spider. Then apply 70% reduction in effect of spider on survival due to hypothesized interference of Met52 with spider. Then apply 50% reduction in survival of remaining ticks; this assumes no interference effect of spider on Met52. For each randomization run, determine P values for effect of Met52, spider, and Met52*spider interaction. As check that these results make sense, also compute delta AIC between model (survival ~ spider + Met52) vs. model (survival ~ spider + Met52 + spider*Met52). Deploy number of microcosms with sample size with &gt;80% power to detect met52*spider interaction effect.

install packages
================

``` r
list.of.packages <- c("lme4", "MuMIn", "dplyr", "ggplot2", "AICcmodavg", "grid", "cowplot", "lubridate")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

print(new.packages)
```

    ## character(0)

``` r
library(lme4)
```

    ## Loading required package: Matrix

``` r
library(MuMIn)#use for getting p values
library("ggplot2")
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library("AICcmodavg")#Use this due to small sample size
```

    ## 
    ## Attaching package: 'AICcmodavg'

    ## The following objects are masked from 'package:MuMIn':
    ## 
    ##     AICc, DIC, importance

    ## The following object is masked from 'package:lme4':
    ## 
    ##     checkConv

``` r
library(grid)
require(cowplot)
```

    ## Loading required package: cowplot

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

``` r
library("lubridate")
```

    ## 
    ## Attaching package: 'lubridate'

    ## The following object is masked from 'package:base':
    ## 
    ##     date

``` r
# Start the clock! (to find out how long this takes)
old <- Sys.time() # get start time

#read in data from previous experiment
D = read.csv("raw_data_files/Tick Predator Project_2015_Tick Survival Data.csv")
#subset for data of interest: with leaf litter in microcosms 
D = subset(D, Litter == "LitRef")

#Assume Met52 causes a 50% reduction in tick survival. 
Met52.reduce.ticks = 0.5#

#Assume Met52 causes a 70% reduction in the effects of spider predation on tick survival. 
Met52.reduce.pred = 0.7

#Create vector of potential sample sizes to examine in power analysis. Set 40 as upper limit; anything more than that may be infeasible logistically. 
sample_size = seq(from = 2, to = 40,by = 2)

#Set number of repetitions in bootstrap power analysis. Initially use a low number e.g. 10 to make sure code works, then increase to 1000. 
reps = 10
out = NULL#Initialize data.frame to hold results from randomization. 
for (a in 1:reps){#for each rep
  #print(a)#print a to know how far along randomization is
  for (b in 1:length(sample_size)){#for each sample size
  N = sample_size[b]#define N with sample size of this randomization run 

  D.PredRef = subset(D, Predator =="PredRef")#subset for no-predator microcosms
  D.PredRef.ticks= D.PredRef$Ixodes#number of ticks that survived without spider
  D.PredAdd = subset(D, Predator =="Addition")  #get tick data for samples with spider
  D.PredAdd.ticks= D.PredAdd$Ixodes#ticks that survived with spider additoin
  Predator = c(rep("PredRef", N*2), rep("PredAdd", N*2) )#label to be used with randomly generated samples with predator
  pred_num = c(rep(1, N*2), rep(2, N*2))#Numeric label for predator status of microcosms
  Pesticide = rep(c(rep("NoPesticide",N), rep("Met52", N)),2)#Initialize vector of labels for Pesticide/No Pesticide
                
  pesticide_num = c(rep(1, N), rep(2, N),
                    rep(1, N), rep(2, N))#Initialize vector of numeric pesticide labels, to be used in modeling (because sometimes weird things happen when words are used)
  tick.rand = rep(NA, N*4)#initialize vector of tick survival values
  
  i.nopred.nopesticide.ref = which(Predator == "PredRef" & Pesticide == "NoPesticide")#Find indices with no predator and no pesticide 
  #To generate tick survival for no spider / no pesticide microcosms, sample with replacement from tick survival in no-spider microcosms.
  tick.rand[i.nopred.nopesticide.ref] = sample(D.PredRef.ticks, size = N, replace = TRUE)
  
  #no predator, yes Met52 
  i.PredRef.Met52 = which(Predator == "PredRef" & Pesticide == "Met52")  #get indices of no spider / Met52 microcosms
  #generate randomized data by sampling with replacement from no-spider microcosms, then apply reduction in survival by multiplying each value by Met52.reduce.ticks
  tick.rand[i.PredRef.Met52] = Met52.reduce.ticks*
    sample(D.PredRef.ticks, size = N, replace = TRUE)

  #predator, no Met52
  #get indices
  i.PredAdd.NoPesticide = which(Predator == "PredAdd" & Pesticide == "NoPesticide")
  #sample with replacement from microcosms with spider
  tick.rand[i.PredAdd.NoPesticide] = sample(D.PredAdd.ticks, size = N, 
                                            replace = TRUE)
  
  #yes predator, yes Met52
  #get a random set of data from PredAdd
  PredAdd.sample = sample(D.PredAdd.ticks, size = N, replace = TRUE)
  #get a random set of data from PredRef 
  PredRef.sample = sample(D.PredRef.ticks, size = N, replace = TRUE)
  #find the difference between PredRef.sample and PredAdd.sample. Assume this is the number killed by predators, relative to number that would have survived in no-spider case
  Pred.effect.rand = PredRef.sample - PredAdd.sample
  i.minus = which(Pred.effect.rand<0)#find indices that have more than 100% mortality
  Pred.effect.rand[i.minus ] = 0#assume there can only be 100% mortality, no more. 
  #Now multiply the reduction in ticks expected w/ predators by a fraction 
  #(Met52.reduce.pred) that reduces the effect of predators on ticks, due to interference of Met52 with spider. 
  Pred.Met52.diffs = Met52.reduce.pred * Pred.effect.rand
  #now find the number of ticks that remain for Met52 to have an effect on, after natural variation and predators have their (reduced by Met52) effect. 
  Pred.remain.for.Met52 = PredRef.sample - Pred.Met52.diffs
  i.PredAdd.Met52 = which(Predator == "PredAdd" & Pesticide == "Met52") 
  tick.rand[i.PredAdd.Met52] = Met52.reduce.ticks*Pred.remain.for.Met52
  D.rand = data.frame(Predator, pred_num, Pesticide, pesticide_num,
                      tick.rand)
  #assign pred_num and pesticide_num as factor
  D.rand$pred_num = factor(D.rand$pred_num)
  D.rand$pesticide_num = factor(D.rand$pesticide_num)
  #fit model that includes effect of pesticide and spider, no interaction
  fit.PP <- glm(tick.rand ~ pred_num + pesticide_num, data =D.rand)
  #fit model that includes effect of predator*pesticide interaction
  fit.PPI <- glm(tick.rand ~ pred_num + pesticide_num + 
                   pred_num * pesticide_num, data =D.rand)
  delta.AIC = (fit.PP$aic- fit.PPI$aic) #get difference in AIC values
  PP.AIC = fit.PP$aic
  PPI.AIC = fit.PPI$aic
  an = summary( fit.PPI)#use analysis of variance to get p values 
  pred.p= an$coefficients[2,4]
  pest.p = an$coefficients[3,4]
  pred.pesticide.int.p = an$coefficients[4,4]
  samples = N
  #add results of this randomization run to temporary data.frame  
  df.tmp = data.frame(samples, pred.p, pest.p, pred.pesticide.int.p,
                      PP.AIC,
                      PPI.AIC,
                      delta.AIC)
  out = rbind(out, df.tmp)
  }
}
#Now get summary of results of randomization
df = out
rep_chk= rep(NA, length(sample_size))
pred.sig05 = rep(NA, length(sample_size))#Initialize vectors for number and fraction of randomization runs significant at P<0.05 level for spider
pred.frac_sig05 = rep(NA, length(sample_size))

#Initialize vector for number, fraction of randomization runs significant at P<0.05 level for pesticide
pest.sig05 = rep(NA, length(sample_size))
pest.frac_sig05 = rep(NA, length(sample_size))

interact.sig05 = rep(NA, length(sample_size))#Initialize vector for number of randomization runs significant at P<0.05 level for spider*Met52 interaction
interact.frac_sig05 = rep(NA, length(sample_size))#Initialize vector for fraction of randomization runs significant
#initialize vector of number, fraction of randomization runsd with delta AIC at least 2
delta.AIC.2plus.sum = rep(NA, length(sample_size))
delta.AIC.2plus.frac = rep(NA, length(sample_size))
for (a in 1:length(sample_size )){#for each sample size
  #number of runs for this sample size
  rep_chk[a] = length(which(df$samples==sample_size[a]))
  interact.sig05[a] = length(which(df$samples==sample_size[a] & 
                                     df$pred.pesticide.int.p<0.05))
  interact.frac_sig05[a] =  interact.sig05[a]/rep_chk[a]
  
  
  pred.sig05[a] = length(which(df$samples==sample_size[a] & 
                            df$pred.p<0.05))
  pred.frac_sig05[a] =  pred.sig05[a]/rep_chk[a]

  pest.sig05[a] = length(which(df$samples==sample_size[a] & 
                                 df$pest.p<0.05))
  pest.frac_sig05[a] =  pest.sig05[a]/rep_chk[a]
  delta.AIC.2plus.sum[a] = length(which(df$samples==sample_size[a] & 
                                          df$delta.AIC>=2))
  delta.AIC.2plus.frac[a] =delta.AIC.2plus.sum[a]/rep_chk[a]
}
df.sum = data.frame(sample_size, 
                    pred.sig05,
                    pred.frac_sig05,
                    pest.sig05,
                    pest.frac_sig05,
                    interact.sig05, 
                    interact.frac_sig05, 
                    delta.AIC.2plus.sum,
                    delta.AIC.2plus.frac,
                    rep_chk )
boot = df.sum
save(boot, file = "Rdata_files/boot.Rdata")
#write summary to file
filename_var = "output/Met52.reduce.pred"
write.csv(df.sum, file = paste(filename_var, Met52.reduce.pred,"reps",reps, 
                               "power_detect.csv", sep = "."),
          row.names = FALSE)
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format
```

    ## Time difference of 0.978723 secs

##### make plot with output of bootstrap

``` r
load("Rdata_files/boot.Rdata")

plot.tmp = ggplot()+
  geom_point(data = boot,
             mapping = aes(x = sample_size,
                           y =interact.frac_sig05))+
                 labs(x = "Replicate of sample size", 
        y="Fraction of randomization runs w/ significant Met52*spider interaction")+
     ylim(0,1)


ggsave(plot.tmp,filename=paste("output/Figure.boot",".jpeg",sep=""),
        dpi = 600,
       width = 6,
       height = 7,
       units = "in")
```

##### How much Met52 should be used in the experiment? What dosage should be applied to the soil cores and to the inside of the organza bag that goes around the core?

###### Decision: Use 3X recommended dosage, based on this being equivalent on per-area basis to dosage used in other papers.

``` r
#example paper: Stafford and Allan: 2.5 X 10^11 spores (conidia) / 100 m^2

#Confirm that 3X dosage is in same order of magnitude as in Stafford and Allan. 
#(3oz/1000 sq ft)*(30g/1oz)*(5*10^10spores/1g)*(1sqft/0.09sqm) 
bag =(3/1000)*(30/1)*(5*10^10/1)*(1/0.092903)
format(bag,scientific = TRUE)
```

    ## [1] "4.843762e+10"

``` r
#(9oz/1000 sq ft)*(30g/1oz)*(5*10^10spores/1g)*(1sqft/0.09sqm) 
core =(9/1000)*(30/1)*(5*10^10/1)*(1/0.092903)
format(core,scientific = TRUE)
```

    ## [1] "1.453129e+11"

##### Analysis of results

###### Read in data on ticks recovered from microcosms and clean it up.

###### Decision: It wasn't possible to collect only gravid or only non-gravid S. ocreata females, so a mix was used. Can data from microcosms with gravid versus non-gravid spider be treated equally in analysis, or is it necessary to account for spider reproductive status in analysis?

###### Approach: determine if there was an effect of spider reproductive status. If no effect of spider reproductive status, then do not include reproductive status in further analysis.

###### Compare fit of model that includes only effect of spider presence vs. model that also includes spider reproductive status.

``` r
load("Rdata_files/C.Rdata")
t1.flat.spid.only <- lm(frac.flat.nymphs.recovered ~ treatment.spider, data = C)


t1.eggs  <- lm(frac.flat.nymphs.recovered ~ 
                     spider.eggs.1, data = C)

x <- c(AICc(t1.flat.spid.only),
        AICc(t1.eggs))

n.param = c(length(t1.flat.spid.only$coefficients),
            length(t1.eggs$coefficients))
#,
 delta <- x - min(x)               # AIC differences
 L <- exp(-0.5 * delta)            # relative likelihoods of models
 model.names <- c(
   "spider",
    "gravid")#

 #residuals
df.residuals <- c(t1.flat.spid.only$df.residual,
       t1.eggs$df.residual)#,

#AIC weights
 w <- L/sum(L)  
 #sort everything by weights
 sorted = sort(w, index.return = TRUE,
               decreasing = TRUE)
 model.names.sorted = model.names[sorted$ix]
 n.param = n.param[sorted$ix]
 AIC.sorted = x[sorted$ix]
 delta.sorted = delta[sorted$ix]
 df.residual.sorted = df.residuals[sorted$ix]
 AIC.weights.sorted = sorted$x
 L = L[sorted$ix]
# 
 df.aic = data.frame(model.names.sorted,
                df.residual.sorted,
                n.param,
                round(AIC.sorted, digits = 2),
                round(delta.sorted, digits = 2),
                round(L, digits = 2),
                round(AIC.weights.sorted, digits = 2)
                )
#make nicer names for table
names(df.aic) = c("model",
              "residual df",
              "number parameters",
              "AICc",
              "delta AIC",
              "Likelihood",
              "AIC weight")

#output table for appendix for manuscript
 write.csv(df.aic, file = "output/Appendix.S1.csv",
           row.names = FALSE)
```

###### There is similar support for model that includes spider presence/absence vs. model that includes reproductive status. Exclude reproductive status from further analysis.

##### What are the effects of Met52, spider, and Met52\*spider interaction on ticks? Decision: Compare fits of alternative models including effects of spider, Met52, or interaction, as predictors of flat and engorged tick survival. Analyze effects of treatment on fraction of flat ticks recovered. Output data for data repository to accompany paper.

``` r
#note need to use glm rather than lm to get aic value as output
load("Rdata_files/C.Rdata")

#Define alternate models, including effects of spider only, Met52 only, both spider and Met52, and spider*Met52 interaction, and null (intercept) model
t1.flat.spid.only <- lm(frac.flat.nymphs.recovered ~ treatment.spider, data = C)

t1.flat.met52.only <- lm(frac.flat.nymphs.recovered ~ treatment.met52, data = C)

t1.flat.met52.spid <- lm(frac.flat.nymphs.recovered ~ treatment.met52 +treatment.spider, 
                             data = C)

t1.flat.met52.spid.int <- lm(frac.flat.nymphs.recovered ~ treatment.met52 + treatment.spider +
                        treatment.met52*treatment.spider, data = C)

t1.null <- lm(frac.flat.nymphs.recovered ~ 1, data = C)

#get AIC values, delta AIC, weights; would be faster to use aictab to do this, but discovered that after-the-fact
 x <- c(
   AICc(t1.flat.spid.only),
        AICc(t1.flat.met52.only),
        AICc(t1.flat.met52.spid),
      AICc(t1.flat.met52.spid.int),
      AICc(t1.null))

 np <- c(
   length(t1.flat.spid.only$coefficients),
        length(t1.flat.met52.only$coefficients),
        length(t1.flat.met52.spid$coefficients),
      length(t1.flat.met52.spid.int$coefficients),
      length(t1.null$coefficients))
       
 delta <- x - min(x)               # AIC differences
 L <- exp(-0.5 * delta)            # relative likelihoods of models
 model.names <- c(
   "spider",
       "Met52",
        "spider + Met52",
        "spider + Met52 + spider*Met52",
       "intercept")
 
df.residuals <- c(t1.flat.spid.only$df.residual,
        t1.flat.met52.only$df.residual,
        t1.flat.met52.spid$df.residual,
      t1.flat.met52.spid.int$df.residual,
      t1.null$df.residual)
       
 w <- L/sum(L)  
 sorted = sort(w, index.return = TRUE,
               decreasing = TRUE)
 model.names.sorted = model.names[sorted$ix]
 np = np[sorted$ix]
 AIC.sorted = x[sorted$ix]
 delta.sorted = delta[sorted$ix]
 df.residual.sorted = df.residuals[sorted$ix]
 AIC.weights.sorted = sorted$x
 L = L[sorted$ix]
# 
 df.aic = data.frame(model.names.sorted,
                df.residual.sorted,
                np,
                round(AIC.sorted, digits = 2),
                round(delta.sorted, digits = 2),
                round(L, digits = 2),
                round(AIC.weights.sorted, digits = 2)
                )
 
names(df.aic) = c("Model",
              "Residual df",
              "Number parameters",
              "AICc",
              "Delta AIC",
              "Likelihood",
              "AIC weight")
 write.csv(df.aic, file = "output/Table1.csv", row.names = FALSE)
 
 #workaround for making table for manuscript for interaction model -- Table 2
tmp = summary(t1.flat.met52.spid.int)
tmp = tmp$coefficients
names(tmp) = c("Coefficient estimate",
               "Coefficient std. error",
               "t value",
               "P(>|t|)")
write.csv(tmp, file = "output/Table.2.csv")
tmp =read.csv("output/Table.2.csv")
names(tmp) = c("Term",
  "Coefficient estimate",
               "Coefficient std. error",
               "t value",
               "P(>|t|)")
tmp$Term = c("intercept",
  "Met52",
       "spider",
       "Met52*spider")
tmp[,c(2:4)] = round(tmp[,c(2:4)], digits = 2)
write.csv(tmp, file = "output/Table.2.csv",
          row.names = FALSE)

#Now output data for repository
keep = c("core.id",
         "treatment.met52",
         "treatment.spider",
  "frac.flat.nymphs.recovered",
         "frac.engorged.nymphs.recovered",
  "live.spider",
  "spider.eggs.1")
C.archive = C[,keep]
names(C.archive) = c("microcosm.id",
         "treatment.met52",
         "treatment.spider",
  "fraction.flat.nymphs.recovered",
         "fraction.engorged.nymphs.recovered",
  "spider.alive.end",
  "spider.reproductive.status.begin")
write.csv(C.archive, file = "output/DataS1.csv",
          row.names = FALSE)
names(C)
```

    ##  [1] "core.id"                                  
    ##  [2] "experiment.title"                         
    ##  [3] "treatment.met52"                          
    ##  [4] "treat.met52.title"                        
    ##  [5] "treat.spider.title"                       
    ##  [6] "treatment.spider"                         
    ##  [7] "bag.check"                                
    ##  [8] "flag.color"                               
    ##  [9] "flag.labeled"                             
    ## [10] "date.core.in"                             
    ## [11] "pvc.labeled"                              
    ## [12] "zone"                                     
    ## [13] "row"                                      
    ## [14] "date.sprayed"                             
    ## [15] "offset.days"                              
    ## [16] "core.removed.planned"                     
    ## [17] "core.removed.date"                        
    ## [18] "core.sorted.date"                         
    ## [19] "time.tot"                                 
    ## [20] "time.sprayed"                             
    ## [21] "bag.sprayed"                              
    ## [22] "flat.nymphs"                              
    ## [23] "flat.date"                                
    ## [24] "fed.nymphs"                               
    ## [25] "fed.date"                                 
    ## [26] "fed.id"                                   
    ## [27] "spider"                                   
    ## [28] "spider.date"                              
    ## [29] "spider.eggs.1"                            
    ## [30] "core.sprayed"                             
    ## [31] "twist.tie"                                
    ## [32] "rubber.band"                              
    ## [33] "zip.tie"                                  
    ## [34] "spider.present.7.26"                      
    ## [35] "notes"                                    
    ## [36] "nymphs.flat.recovered"                    
    ## [37] "nymphs.engorged.recovered"                
    ## [38] "larvae.recovered"                         
    ## [39] "live.spider"                              
    ## [40] "non.s.ocreata.spider.0.5.cm"              
    ## [41] "spider.eggs"                              
    ## [42] "spiderlings.present"                      
    ## [43] "non.s.ocreata.1.cm"                       
    ## [44] "spider.juveniles.present"                 
    ## [45] "adult.male.recovered"                     
    ## [46] "adult.female.recovered"                   
    ## [47] "berlese.flat.nymphs.recovered"            
    ## [48] "cent.0.5.cm.over"                         
    ## [49] "cent.over.1.cm"                           
    ## [50] "beetle.over.0.5.cm"                       
    ## [51] "beetle.over.1.cm"                         
    ## [52] "berlese.adult.s.ocreata.recovered"        
    ## [53] "berlese.juve.s.ocreata"                   
    ## [54] "berlese.adult.scapularis.recovered"       
    ## [55] "berlese.adult.male.scapularis.recovered"  
    ## [56] "berlese.adult.female.scapularis.recovered"
    ## [57] "live.spider.numeric"                      
    ## [58] "frac.flat.nymphs.recovered"               
    ## [59] "frac.engorged.nymphs.recovered"           
    ## [60] "treatment.met52.char"                     
    ## [61] "treatment.spider.char"                    
    ## [62] "treatment.char"                           
    ## [63] "gravid.char"                              
    ## [64] "gravid.met52.char"                        
    ## [65] "gravid.binary"

###### analyze effects of treatment on engorged tick recovery

``` r
load("Rdata_files/C.Rdata")

C$y = C$frac.engorged.nymphs.recovered
C = subset(C, !is.na(frac.engorged.nymphs.recovered))

#make alternative models
t1.eng.spid.only <- lm(y ~ treatment.spider, data = C)
summary(t1.eng.spid.only)
```

    ## 
    ## Call:
    ## lm(formula = y ~ treatment.spider, data = C)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.2273 -0.2273 -0.1744  0.2727  1.2727 
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        0.22727    0.04951    4.59 1.52e-05 ***
    ## treatment.spider1 -0.05285    0.07043   -0.75    0.455    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3284 on 85 degrees of freedom
    ## Multiple R-squared:  0.006582,   Adjusted R-squared:  -0.005105 
    ## F-statistic: 0.5632 on 1 and 85 DF,  p-value: 0.4551

``` r
t1.eng.met52.only <- lm(y ~ treatment.met52, data = C)
summary(t1.eng.met52.only)
```

    ## 
    ## Call:
    ## lm(formula = y ~ treatment.met52, data = C)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.2234 -0.2234 -0.1750  0.2766  1.2766 
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)       0.17500    0.05196   3.368  0.00114 **
    ## treatment.met521  0.04840    0.07069   0.685  0.49539   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3286 on 85 degrees of freedom
    ## Multiple R-squared:  0.005485,   Adjusted R-squared:  -0.006215 
    ## F-statistic: 0.4688 on 1 and 85 DF,  p-value: 0.4954

``` r
t1.eng.met52.spid <- lm(y ~ treatment.met52 +treatment.spider, 
                             data = C)

t1.eng.met52.spid.int <- lm(y ~ treatment.met52 + treatment.spider +
                        treatment.met52*treatment.spider, data = C)
summary(t1.eng.met52.spid.int)
```

    ## 
    ## Call:
    ## lm(formula = y ~ treatment.met52 + treatment.spider + treatment.met52 * 
    ##     treatment.spider, data = C)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.2292 -0.2250 -0.1250  0.2750  1.2708 
    ## 
    ## Coefficients:
    ##                                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                         0.225000   0.073949   3.043  0.00314
    ## treatment.met521                    0.004167   0.100127   0.042  0.96691
    ## treatment.spider1                  -0.100000   0.104580  -0.956  0.34174
    ## treatment.met521:treatment.spider1  0.088225   0.142299   0.620  0.53696
    ##                                      
    ## (Intercept)                        **
    ## treatment.met521                     
    ## treatment.spider1                    
    ## treatment.met521:treatment.spider1   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3307 on 83 degrees of freedom
    ## Multiple R-squared:  0.0165, Adjusted R-squared:  -0.01905 
    ## F-statistic: 0.4641 on 3 and 83 DF,  p-value: 0.7081

``` r
t1.eng.null <- lm(y ~ 1, data = C)
summary(t1.eng.null)
```

    ## 
    ## Call:
    ## lm(formula = y ~ 1, data = C)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.2011 -0.2011 -0.2011  0.2989  1.2989 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.20115    0.03512   5.727 1.47e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3276 on 86 degrees of freedom

``` r
 x <- c(
   AICc(t1.eng.spid.only),
        AICc(t1.eng.met52.only),
        AICc(t1.eng.met52.spid),
      AICc(t1.eng.met52.spid.int),
      AICc(t1.eng.null))
  
 np <- c(
   length(t1.eng.spid.only$coefficients),
        length(t1.eng.met52.only$coefficients),
        length(t1.eng.met52.spid$coefficients),
      length(t1.eng.met52.spid.int$coefficients),
      length(t1.eng.null$coefficients))
    
 delta <- x - min(x)               # AIC differences
 L <- exp(-0.5 * delta)            # relative likelihoods of models
 model.names <- c(
   "spider",
       "Met52",
        "spider + Met52",
        "spider + Met52 + spider*Met52",
       "intercept")
 
df.residuals <- c(t1.eng.spid.only$df.residual,
        t1.eng.met52.only$df.residual,
        t1.eng.met52.spid$df.residual,
      t1.eng.met52.spid.int$df.residual,
      t1.eng.null$df.residual)
  
 w <- L/sum(L)  
 sorted = sort(w, index.return = TRUE,
               decreasing = TRUE)
 model.names.sorted = model.names[sorted$ix]
 np = np[sorted$ix]
 AIC.sorted = x[sorted$ix]
 delta.sorted = delta[sorted$ix]
 df.residual.sorted = df.residuals[sorted$ix]
 AIC.weights.sorted = sorted$x
 L = L[sorted$ix]
# 
 df.aic = data.frame(model.names.sorted,
                df.residual.sorted,
                np,
                round(AIC.sorted, digits = 2),
                round(delta.sorted, digits = 2),
                round(L, digits = 2),
                round(AIC.weights.sorted, digits = 2)
                )
 
names(df.aic) = c("Model",
              "Residual df",
              "Number parameters",
              "AICc",
              "Delta AIC",
              "Likelihood",
              "AIC weight")
 write.csv(df.aic, file = "output/Table.3.eng.spider.Met52.AIC.csv", row.names = FALSE)
```

###### graph: boxplots of proportion of flat tick and engorged tick survival by treatment

    ## Warning: Removed 2 rows containing non-finite values (stat_boxplot).

##### Is spider survival reduced by Met52 treatment?

###### Decision: Compare model with Met52 vs. null (intercept) model

``` r
load("Rdata_files/C.Rdata")

C.spider = subset(C, treatment.spider==1)
t1.met52 <- glm(live.spider ~ treatment.met52, data = C.spider,
                 family = "binomial")
summary(t1.met52)
```

    ## 
    ## Call:
    ## glm(formula = live.spider ~ treatment.met52, family = "binomial", 
    ##     data = C.spider)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.5518  -0.7775  -0.7775   0.8446   1.6394  
    ## 
    ## Coefficients:
    ##                  Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)        0.8473     0.4879   1.736  0.08249 . 
    ## treatment.met521  -1.8888     0.6809  -2.774  0.00554 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 59.401  on 42  degrees of freedom
    ## Residual deviance: 50.837  on 41  degrees of freedom
    ## AIC: 54.837
    ## 
    ## Number of Fisher Scoring iterations: 4

``` r
t1.null = glm(live.spider ~ 1, data = C.spider,
                 family = "binomial")
 x <- c(
   AICc(t1.met52),
   AICc(t1.null))
 
 np <- c(
   length(t1.met52$coefficients),
        length(t1.null$coefficients))

      
 delta <- x - min(x)               # AIC differences
 L <- exp(-0.5 * delta)            # relative likelihoods of models
 model.names <- c(
   "Met52",
       "intercept")
 
 df.residuals <- c(t1.met52$df.residual,
       t1.null$df.residual)
       
 w <- L/sum(L)  
 sorted = sort(w, index.return = TRUE,
               decreasing = TRUE)
 model.names.sorted = model.names[sorted$ix]
 AIC.sorted = x[sorted$ix]
 delta.sorted = delta[sorted$ix]
 np = np[sorted$ix]
 
 df.residual.sorted = df.residuals[sorted$ix]
 AIC.weights.sorted = sorted$x
 L = L[sorted$ix]
# 
 df.aic = data.frame(model.names.sorted,
                df.residual.sorted,
                np,
                round(AIC.sorted, digits = 2),
                round(delta.sorted, digits = 2),
                round(L, digits = 2),
                round(AIC.weights.sorted, digits = 2)
                )
# 
names(df.aic) = c("model",
              "residual df",
              "number parameters",
              "AICc",
              "delta AIC",
              "Likelihood",
              "AIC weight")

 write.csv(df.aic, file = "output/Table4.spider.survival.Met52.AIC.csv", row.names = FALSE)
```

###### graph spider survival by treatment

``` r
load("Rdata_files/C.Rdata")

C = subset(C, treatment.spider == 1)
fs= 18
xlabs = c("H2O", "Met52")
spider.plot = ggplot()+
    geom_boxplot(data=C, 
             mapping=aes(x=treatment.char, 
                         y = live.spider.numeric))+
  labs(x = "treatment", 
       y="spider survival")+
    ylim(0,1)+
    scale_x_discrete(labels = xlabs)+
     theme(axis.text=element_text(size=fs))+
  theme(axis.title = element_text(size = fs))
save_plot( "output/Figure.spider.survival.jpg", spider.plot, 
           dpi = 600)
```

##### Do spiders or Met52 affect tick questing?

###### Decision: Use model comparison to test for effects of treatment and spider activity on number of ticks questing for hosts (when I stick my face next to the microcosm).

###### Clean up questing data and output for repository

``` r
load("Rdata_files/C.Rdata")

#read in data on ticks questing
Q = read.csv(file = "raw_data_files/tritrophic cores - questing_DATA.csv")
#limit to the cores we are using in survival analysis (remove cores that I made mistake with spray application with)
names(Q)= tolower(names(Q))
Q = subset(Q, core.id %in% C$core.id)
#get the number of questing assessments done 
names(Q)[names(Q)=="live.spider"]="spider.alive.questing"

#merge questing data with rest of data
O = merge(C,Q, by = "core.id")
O$spider.alive.questing.numeric = O$spider.alive.questing
O$spider.alive.questing = factor(O$spider.alive.questing)

#fix date
O$date = paste(O$date,"/2017", sep = "")
O$date <- as.Date(O$date, "%m/%d/%Y")
out = NULL

#for each core get the most recent questing assessment
core.ids = unique(O$core.id)
for (a in 1:length(core.ids)){
 tmp = subset(O, core.id== core.ids[a])
 tmp = subset(tmp, date == max(tmp$date))
 out = rbind(out, tmp)
}
O = out

# treatment.met52.char = rep(NA, dim(O)[1])
# treatment.met52.char[O$treatment.met52==0]="H2O"
# treatment.met52.char[O$treatment.met52==1]="Met52"
# O$treatment.met52.char = treatment.met52.char
# 
# treatment.spider.char = rep(NA, dim(O)[1])
# treatment.spider.char[O$treatment.spider==0]="no spider"
# treatment.spider.char[O$treatment.spider==1]="spider"
# O$treatment.spider.char = treatment.spider.char
# 
# O$treatment.char = paste(O$treatment.met52.char, 
#                          O$treatment.spider.char,
#                          sep = "_")
# O = subset(O, treatment.char !="NA_NA")

#get mean number of ticks questing across the up-to-four 30-second samples
ticks = rep(NA, dim(O)[1])
for (a in 1:length(ticks)){
  ticks[a] = mean(O$x1[a],
              O$x2[a],
              O$x3[a],
              O$x4[a],
              O$x5[a],
              na.rm = TRUE)
}
ticks[ticks==-Inf]= NA
O$ticks = ticks
O$frac.nymphs.questing = O$ticks/O$nymphs.flat.recovered
O$frac.nymphs.questing[O$frac.nymphs.questing == Inf] = NA
O$frac.nymphs.questing[O$frac.nymphs.questing == "NaN"] = NA

save(O, file = "Rdata_files/O.Rdata")

keep = c("core.id",
         "frac.nymphs.questing",
         "treatment.met52",
         "treatment.spider",
         "spider.alive.questing",
         "live.spider")

#output data for repository. 
O.archive = O[,keep]
names(O.archive) =  c("microcosm.id",
         "fraction.nymphs.questing",
         "treatment.met52",
         "treatment.spider",
         "spider.observed",
         "spider.alive.end")
write.csv(O.archive, file = "output/Data.S2.csv")
#output for figshare
```

##### run models of tick questing re: Met52, spider activity, spider survival

``` r
load("Rdata_files/O.Rdata")

t1.met52 = lm(frac.nymphs.questing ~ treatment.met52, data = O)

t1.spider <- lm(frac.nymphs.questing ~ treatment.spider, data = O)

t1.active.spider <- lm(frac.nymphs.questing ~ spider.alive.questing, data = O)

t1.spider.live <- lm(frac.nymphs.questing ~ 
                                live.spider, data = O)

t1.met52.spid = lm(frac.nymphs.questing ~
                     treatment.met52+
                     treatment.spider, 
                   data = O)

t1.met52.live = lm(frac.nymphs.questing ~
                     treatment.met52+
                     live.spider, 
                   data = O)

t1.met52.quest = lm(frac.nymphs.questing ~
                     treatment.met52+
                     spider.alive.questing, 
                   data = O)

t1.spider.live <- lm(frac.nymphs.questing ~ 
                       treatment.spider + live.spider,
                     data = O)

t1.spider.quest <- lm(frac.nymphs.questing ~
                        treatment.spider + spider.alive.questing,
                      data = O)

t1.spider.live.quest <- lm(frac.nymphs.questing ~ 
                                live.spider + spider.alive.questing, data = O)

t1.all.spider <- lm(frac.nymphs.questing ~ 
                                treatment.spider+
                      live.spider + 
                      spider.alive.questing,
                    data = O)

t1.spider.quest.met52 <- lm(frac.nymphs.questing ~ 
                              treatment.spider+  
                              spider.alive.questing + 
                                  treatment.met52, data = O)


t1.spider.live.met52 <- lm(frac.nymphs.questing ~ 
                              treatment.spider+  
                              live.spider + 
                                  treatment.met52, data = O)

t1.spider.live.quest.met52 <- lm(frac.nymphs.questing ~ 
                                treatment.spider +
                                   live.spider + 
                                  spider.alive.questing+
                                  treatment.met52, data = O)

t1.null <- lm(frac.nymphs.questing ~ 1, data = O)

#this could be done faster with model.sel
 x <- c(AICc(t1.met52),
   AICc(t1.spider),
        AICc(t1.active.spider),
        AICc(t1.spider.live),
      AICc(t1.met52.spid ),
      AICc(t1.met52.live ),
      AICc(t1.met52.quest ),
      AICc(t1.spider.live ),
      AICc(t1.spider.quest ),
      AICc(t1.spider.live.quest ),
      AICc(t1.all.spider ),
      AICc(t1.spider.quest.met52 ),
      AICc(t1.spider.live.met52 ),
      AICc(t1.spider.live.quest.met52 ),
      AICc(t1.null))
 
  np <- c(length(t1.met52$coefficients),
   length(t1.spider$coefficients),
        length(t1.active.spider$coefficients),
        length(t1.spider.live$coefficients),
      length(t1.met52.spid$coefficients ),
      length(t1.met52.live$coefficients ),
      length(t1.met52.quest$coefficients ),
      length(t1.spider.live$coefficients ),
      length(t1.spider.quest$coefficients ),
      length(t1.spider.live.quest$coefficients ),
      length(t1.all.spider$coefficients ),
      length(t1.spider.quest.met52$coefficients ),
      length(t1.spider.live.met52$coefficients ),
      length(t1.spider.live.quest.met52$coefficients ),
      length(t1.null$coefficients))
  
  df.residual <- c((t1.met52$df.residual),
   (t1.spider$df.residual),
        (t1.active.spider$df.residual),
        (t1.spider.live$df.residual),
      (t1.met52.spid$df.residual ),
      (t1.met52.live$df.residual ),
      (t1.met52.quest$df.residual ),
      (t1.spider.live$df.residual ),
      (t1.spider.quest$df.residual ),
      (t1.spider.live.quest$df.residual ),
      (t1.all.spider$df.residual ),
      (t1.spider.quest.met52$df.residual ),
      (t1.spider.live.met52$df.residual ),
      (t1.spider.live.quest.met52$df.residual ),
      (t1.null$df.residual))
 
 
delta <- x - min(x)               # AIC differences
 L <- exp(-0.5 * delta)            # relative likelihoods of models
model.names = c("t1.met52",
  "t1.spider",
        "t1.active.spider",
        "t1.spider.live",
      "t1.met52.spid",
      "t1.met52.live",
    "t1.met52.quest",
      "t1.spider.live",
      "t1.spider.quest",
      "t1.spider.live.quest",
      "t1.all.spider",
      "t1.spider.quest.met52",
      "t1.spider.live.met52",
      "t1.spider.live.quest.met52",
      "t1.null")
  
model.names.long = c("Met52",
 "spider treatment",
        "spider active",
        "spider survive",
      "Met52 + spider treatment",
      "Met52 + spider lived",
    "Met52 + spider active",
      "spider treatment + spider lived",
      "spider treatment + spider active",
      "spider active + spider lived",
      "spider treatment + spider active + spider lived",
      "Met52 + spider treatment + spider active",
      "Met52 + spider treatment + spider lived",
      "Met52 + spider treatment + spider active + spider lived",
      "Null (intercept)")

 w <- L/sum(L)  
 sorted = sort(w, index.return = TRUE,
               decreasing = TRUE)
 model.names.sorted = model.names[sorted$ix]
  model.names.sorted.long = model.names.long[sorted$ix]
 AIC.sorted = x[sorted$ix]
 np = np[sorted$ix]
 delta.sorted = delta[sorted$ix]
 df.residual.sorted = df.residual[sorted$ix]
 AIC.weights.sorted = sorted$x
 L = L[sorted$ix]
# 
 df.aic = data.frame(#model.names.sorted,
                     model.names.sorted.long,
                 df.residual.sorted,
                 np,
                round(AIC.sorted, digits = 2),
                round(delta.sorted, digits = 2),
                round(L, digits = 2),
                round(AIC.weights.sorted, digits = 2)
                )
# 
names(df.aic) = c(#"Model short",
                  "Model",
              "Residual df",
              "Number parameters",
              "AICc",
              "Delta AIC",
              "Likelihood",
              "AIC weight")
df.aic
```

    ##                                                      Model Residual df
    ## 1          spider treatment + spider active + spider lived          79
    ## 2                         spider treatment + spider active          80
    ## 3                             spider active + spider lived          80
    ## 4  Met52 + spider treatment + spider active + spider lived          78
    ## 5                 Met52 + spider treatment + spider active          79
    ## 6                                            spider active          81
    ## 7                                         Null (intercept)          82
    ## 8                                         spider treatment          81
    ## 9                                    Met52 + spider active          80
    ## 10                                                   Met52          81
    ## 11                                Met52 + spider treatment          80
    ## 12                                          spider survive          80
    ## 13                         spider treatment + spider lived          80
    ## 14                                    Met52 + spider lived          80
    ## 15                 Met52 + spider treatment + spider lived          79
    ##    Number parameters  AICc Delta AIC Likelihood AIC weight
    ## 1                  4 88.00      0.00       1.00       0.31
    ## 2                  3 88.92      0.92       0.63       0.20
    ## 3                  3 88.94      0.94       0.62       0.19
    ## 4                  5 90.20      2.20       0.33       0.10
    ## 5                  4 90.55      2.55       0.28       0.09
    ## 6                  2 92.98      4.98       0.08       0.03
    ## 7                  1 93.22      5.22       0.07       0.02
    ## 8                  2 93.64      5.64       0.06       0.02
    ## 9                  3 95.01      7.01       0.03       0.01
    ## 10                 2 95.37      7.37       0.03       0.01
    ## 11                 3 95.84      7.84       0.02       0.01
    ## 12                 3 95.85      7.85       0.02       0.01
    ## 13                 3 95.85      7.85       0.02       0.01
    ## 14                 3 97.04      9.04       0.01       0.00
    ## 15                 4 98.10     10.10       0.01       0.00

``` r
tmp = summary(t1.all.spider)
#term = names(names(coef(t1.all.spider)))

tmp = tmp$coefficients
names(tmp) = c("Coefficient estimate",
               "Coefficient std. error",
               "t value",
               "P(>|t|)")

write.csv(tmp, file = "Rdata_files/Table.6.csv")

tmp =read.csv("Rdata_files/Table.6.csv")

names(tmp) = c("Term",
  "Coefficient estimate",
               "Coefficient std. error",
               "t value",
               "P(>|t|)")

tmp$Term = c("Intercept",
             "spider treatment",
             "spider lived",
             "spider active")
tmp[,c(2:4)] = round(tmp[,c(2:4)], digits = 2)
tmp[,c(5)] = round(tmp[,c(5)], digits = 3)#P
write.csv(tmp, file = "Rdata_files/Table.6.csv",
          row.names = FALSE)

 write.csv(df.aic, file = "output/quest.spider.Met52.AIC.csv", row.names = FALSE)
```

##### make error bar graph for ticks questing re: treatment

``` r
load("Rdata_files/O.Rdata")
O$live.spider.plot = O$live.spider.numeric
O$live.spider.plot[O$treatment.spider==0]=-1

O$spider.alive.questing.plot = O$spider.alive.questing.numeric
O$spider.alive.questing.plot[O$treatment.spider==0]=-1

T.table = aggregate(x = O$frac.nymphs.questing,
          FUN = mean,
          by = list(O$spider.alive.questing.plot,
                    O$live.spider.plot),
                    na.rm = TRUE)
names(T.table) = c("spider.alive.questing",
                   "spider.survived",
                   "fraction.ticks.seen")

T.table.sd =aggregate(x = O$frac.nymphs.questing,
          FUN = sd,
          by = list(O$spider.alive.questing.plot,
                    O$live.spider.plot),
                    na.rm = TRUE)
names(T.table.sd) = c("spider.alive.questing",
                   "spider.survived",
                   "x")

T.table.se = T.table.sd
T.table.se$x = T.table.se$x/sqrt(dim(O)[1])
T.table$upper = T.table$fraction.ticks.seen + T.table.se$x
T.table$lower = T.table$fraction.ticks.seen - T.table.se$x
spider.visible = rep(NA, dim(T.table)[1])
spider.visible[T.table$spider.alive.questing==-1]="no spider"
spider.visible[T.table$spider.alive.questing==1]="seen"
spider.visible[T.table$spider.alive.questing==0]="not seen"

spider.survived = rep(NA, dim(T.table)[1])
spider.survived[T.table$spider.survived==-1]="no spider"
spider.survived[T.table$spider.survived==0]="died"
spider.survived[T.table$spider.survived==1]="lived"

spider.char = paste(spider.visible, spider.survived, sep = " |\n")
ind = which(spider.char == "no spider |\nno spider")
spider.char[ind] = "no spider \ntreatment"
T.table$spider.char =spider.char 

#plot
plots =  ggplot()+
     geom_errorbar(data=T.table, 
    mapping=aes(x=spider.char, 
                        ymin = T.table$lower,
                ymax = T.table$upper))+
  geom_point(data=T.table, 
    mapping=aes(x=spider.char,
                y = fraction.ticks.seen))
#size = ptsize, shape=21, fill="white" )+
  labs(x = "spider seen in microcosm | spider survived experiment", 
        y="ticks questing as fraction of ticks surviving +-SE")
```

    ## $x
    ## [1] "spider seen in microcosm | spider survived experiment"
    ## 
    ## $y
    ## [1] "ticks questing as fraction of ticks surviving +-SE"
    ## 
    ## attr(,"class")
    ## [1] "labels"

``` r
  ylim(0,1.1)
```

    ## <ScaleContinuousPosition>
    ##  Range:  
    ##  Limits:    0 --  1.1

``` r
   # +
  ggsave(plots,filename=paste("output/Figure.questing",".jpeg",sep=""),
        dpi = 600,
       width = 6,
       height = 7,
       units = "in")
```

##### make boxplot graph for ticks questing re: treatment: no spider / not seen + died / not seen + lived / seen

``` r
load("Rdata_files/O.Rdata")
O$live.spider.plot = O$live.spider.numeric
O$live.spider.plot[O$treatment.spider==0]=-1#never spider

spider.status = rep(NA, dim(O)[1])

O$spider.alive.questing.plot = O$spider.alive.questing.numeric#questing =1 = seen; 0 = not seen
O$spider.alive.questing.plot[O$treatment.spider==0]=-1

spider.status[O$live.spider.plot == -1]= "no \n spider \n treatment"
spider.status[O$treatment.spider == 1 & O$spider.alive.questing.plot == 1]= "spider \n active"
spider.status[O$treatment.spider == 1 & O$live.spider.plot == 1 &
                O$spider.alive.questing.plot == 0]= "spider \n inactive"
spider.status[O$treatment.spider == 1 & O$live.spider.plot == 0] = "spider \n dead"

O$spider.status = spider.status

neworder = c("spider \n inactive", "spider \n active", "spider \n dead", "no \n spider \n treatment")
O <- arrange(transform(O,
             spider.status=factor(spider.status,levels=neworder)),spider.status)


require(cowplot)
plot.spider = ggplot()+
 geom_boxplot(data=O, 
             mapping=aes(x=spider.status, 
                         y = frac.nymphs.questing))+
  labs(x = "spider status at observation time", 
       y="proportion nymphs questing")+
    ylim(0,1)

#p = plot_grid(plot.spider, nrow = 1, align = "v")
save_plot( "output/Figure.3.jpg", plot.spider, 
           dpi = 600)
```

    ## Warning: Removed 9 rows containing non-finite values (stat_boxplot).
