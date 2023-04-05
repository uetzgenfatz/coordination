#' ---
#' title: Disjuntive Agreement
#' author: 
#'  - Phillip M. Alday, Oliver Schallert
#' date: 2022-04-24
#' ---

library("groundhog")


pkgs <- c("lme4", "car", "tidyverse", "here", "emmeans", "Hmisc")

groundhog.library(pkgs, "2022-04-01")

dat <- read.csv("/Users/uetzelsche/Desktop/data_grammatik_2022-04-22_15-52.csv", sep="\t", fileEncoding="UTF-16LE") 

# drop extra rows
dat <- dat[dat$SERIAL != "", ]

# drop some unused columns
cols <- c("CASE", "REF", "MODE", "STARTED", "RA01_CP", "RA01", "TIME001", 
          "TIME002", "TIME003", "TIME004", "TIME005", "TIME006", "TIME007", 
          "TIME_SUM", "MAILSENT",  "LASTDATA", "FINISHED", "Q_VIEWER", 
          "LASTPAGE", "MAXPAGE", "MISSING", "MISSREL", "TIME_RSI", "DEG_TIME")
dat <- select(dat, !cols)

library(dplyr)

# rename some columns
dat <- rename(dat, Liste = QUESTNNR, Informant = SERIAL)

# drop demographic questions
dat <- select(dat, !starts_with("DE"))

# drop practice questions
dat <- select(dat, !starts_with("AW"))

# drop filler 
dat <- select(dat, !starts_with("FI"))

# move to long format
# S[A-E]* 
dat <- pivot_longer(dat, starts_with("S"), names_to = "Bedingung", values_to = "Bewertung")

# Latin-square esque design so we have some NAs to get rid of
dat <- drop_na(dat)

# split out the item ("lexical material") from the condition
dat$Item <- substr(dat$Bedingung, 2, 2)

dat$Bedingung <- substr(dat$Bedingung, 3, 4) 
# get the number + letter because different conditions appeared across items
letter2number <- list(A=1, B=2, C=3, D=4, E=5) # there are more efficient ways to do this, but involve too many tricks
item_nr <- unlist(letter2number[dat$Item], use.names = FALSE)
bedingung_nr <- (item_nr - 1) * 12 + as.numeric(dat$Bedingung) # modular arithmetic to the rescue

# these are just the conditions in the order from the variables_*.csv
bedingungen <- c("1-3-SV1", "1-3-SV3", "1-3-SVpl", "3-1-SV1", "3-1-SV3", 
                 "3-1-SVpl", "1-3-V1S", "1-3-V3S", "1-3-VplS", "3-1-V1S", 
                 "3-1-V3S", "3-1-VplS", "1-3-SV1syn", "1-3-SV3syn", "1-3-SVpl.syn", 
                 "3-1-SV1syn", "3-1-SV3syn", "3-1-SVpl.syn", "1-3-V1synS", 
                 "1-3-V3synS", "1-3-Vpl.synS", "3-1-V1synS", "3-1-V3synS", 
                 "3-1-Vpl.synS", "2-3-SV2", "2-3-SV3", "2-3-SVpl", "3-2-SV2", 
                 "3-2-SV3", "3-2-SVpl", "2-3-V2S", "2-3-V3S", "2-3-VplS", 
                 "3-2-V2S", "3-2-V3S", "3-2-VplS", "2pl-3-SV2pl", "2pl-3-SV3", 
                 "2pl-3-SVpl", "3-2pl-SV2pl", "3-2pl-SV3", "3-2pl-SVpl", 
                 "2pl-3-V2plS", "2pl-3-VS3", "2pl-3-VSpl", "3-2pl-V2plS", 
                 "3-2pl-V3S", "3-2pl-VplS", "2pl-3-SV2pl.syn", "2pl-3-SV3syn", 
                 "2pl-3-SVpl.syn", "3-2pl-SV2pl.syn", "3-2pl-SV3syn", 
                 "3-2pl-SVpl.syn", "2pl-3-V2pl.synS", "2pl-3-V3synS", 
                 "2pl-3-Vpl.synS", "3-2pl-V2pl.synS", "3-2pl-V3synS", 
                 "3-2pl-Vpl.synS" )

dat$Bedingung <- bedingungen[bedingung_nr]

# split Bedingung into three columns
dat <- separate(dat, "Bedingung", into=c("Arg1", "Arg2", "Verb_Wortstellung"), sep="-")

dat$Wortstellung <- ifelse(substr(dat$Verb_Wortstellung, 1, 1) == "S", "SV", "VS")
dat$VerbnahPerson <- ifelse(dat$Wortstellung == "SV", dat$Arg2, dat$Arg1) 
dat$VerbfernPerson <- ifelse(dat$Wortstellung == "SV", dat$Arg1, dat$Arg2) 
dat$Verb <- str_replace_all(dat$Verb_Wortstellung, "[VS]", "")

dat$Wortstellung <- factor(dat$Wortstellung)
dat$VerbnahPerson <- factor(dat$VerbnahPerson) 
dat$VerbfernPerson <- factor(dat$VerbfernPerson) 
dat$Verb <- factor(dat$Verb)
dat$Informant <- factor(dat$Informant)
dat$Item <- factor(dat$Item)

# set some contrasts
contrasts(dat$Wortstellung) <- contr.Sum(levels(dat$Wortstellung))
contrasts(dat$VerbnahPerson) <- contr.Sum(levels(dat$VerbnahPerson)) 
contrasts(dat$VerbfernPerson) <- contr.Sum(levels(dat$VerbfernPerson)) 
contrasts(dat$Verb) <- contr.Sum(levels(dat$Verb))

#' this weird nesting syntax (/) sets up the contrasts so that you essentially have submodels
#' for each of the relevant Verb-Agreements.
#' there are no interactions between VerbnahPerson and VerbfernPerson because the design doesn't really allow for it 
#' Wortstellung is allowed to interact with each of them though -- the formula notation follows the distributive rule.
#'
#' no random slopes for item -- too few items to do anything complex, 
#' especially since not all conditions occur in all items
#' 
#' for informant, we use a different parameterization in the random effects than the fixed effects
#' which has some implications for interpreting them in dept, but does a pretty good job of 
#' allowing for all plausible sources of interindividual variation that would be observable in this experiment
#' when we take effect size into account
#' also....
#' We originally had
#' (1 + Verb + Wortstellung + VerbnahPerson + VerbfernPerson | Informant), 
#' but after doing diagnostics with rePCA() below, we simplified a bit
#' this model is still singular in both Informant and Item, but that might change with additional
#' data, so we’re leaving them in for now.
model <- lmer(Bewertung ~ 1 + Verb / (Wortstellung * (VerbnahPerson + VerbfernPerson)) + 
                        (1|Item) + 
                        (1 + Verb + Wortstellung + VerbnahPerson | Informant), 
              dat, REML=FALSE)
summary(rePCA(model))
summary(model)
Anova(model)

#' for some reason, effects really struggles with the rank deficient fix effects
#' but we can still do some plotting of the model fit. it loooks like the model 
#' does a great job
dat$FittedBewertung <- fitted(model)


# Save plots

model_fit <-
ggplot(dat, aes(x=VerbnahPerson, y=Bewertung)) + 
    stat_summary(aes(color="observed"), fun.data=mean_cl_boot) +
    stat_summary(aes(y=FittedBewertung, color="model"), fun.data=mean_cl_boot) + facet_grid(Wortstellung ~ Verb) + xlab("Verbform") + ylab("Bewertung") + scale_color_discrete(name  ="Modellgüte",
                              breaks=c("model", "observed"),
                              labels=c("Modell", "Daten"))    
    
all_effects <- 
    ggplot(dat, aes(x=Verb, y=Bewertung, color=Wortstellung)) + stat_summary(fun.data=mean_cl_boot) + facet_wrap(~VerbnahPerson) + xlab("Verbform") + ylab("Bewertung") + ylim(0, 4)

single_effects_1 <- 
ggplot(subset(dat, VerbnahPerson=='1'), aes(x=Verb, y=Bewertung, color=Wortstellung)) + stat_summary(fun.data=mean_cl_boot) + xlab("Verbform") + ylab("Bewertung") + ylim(0, 4)

single_effects_2 <- 
    ggplot(subset(dat, VerbnahPerson=='2'), aes(x=Verb, y=Bewertung, color=Wortstellung)) + stat_summary(fun.data=mean_cl_boot) + xlab("Verbform") + ylab("Bewertung") + ylim(0, 4)

single_effects_3 <- 
    ggplot(subset(dat, VerbnahPerson=='3'), aes(x=Verb, y=Bewertung, color=Wortstellung)) + stat_summary(fun.data=mean_cl_boot) + xlab("Verbform") + ylab("Bewertung") + ylim(0, 4)

single_effects_2pl <- 
    ggplot(subset(dat, VerbnahPerson=='2pl'), aes(x=Verb, y=Bewertung, color=Wortstellung)) + stat_summary(fun.data=mean_cl_boot) + xlab("Verbform") + ylab("Bewertung") + ylim(0, 4)

setwd("/Users/uetzelsche/Desktop")

single_effects_1
quartz.save("1sg.png","png", dpi = 150) 

single_effects_2
quartz.save("2sg.png","png", dpi = 150)

single_effects_3
quartz.save("3sg.png","png", dpi = 150) 

single_effects_2pl
quartz.save("2pl.png","png", dpi = 150) 

model_fit
quartz.save("modelfit.png","png", dpi = 150) 

all_effects
quartz.save("alleffects.png","png", dpi = 150) 

