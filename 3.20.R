
library(lmtest)
library(DHARMa)
library(lme4)
library(geoR)
library(spacetime)
library(gstat)
library(sf)
library(sp)
library(nlme)
library(MASS)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gdm)
library(vegan)
library(glmmTMB)


Omastersheet<-read.csv("OriginalMastersheet3.19.csv", header=TRUE)

colnames(Omastersheet)
mastersheet<-Omastersheet[,-c(1,4,8:10,13:21,24:209)] #removing unnecessary columns

#grouping parasites of the same species from diff locations into a single column
clav.mariae<-mastersheet$Clavinema_mariae.gill+mastersheet$Clavinema_mariae.kidney+
  mastersheet$Clavinema_mariae.gonad+
  mastersheet$Clavinema_mariae.dorsalfin+mastersheet$Clavinema_mariae.analfin+
  mastersheet$Clavinema_mariae.cyst.dorsal+
  mastersheet$Clavinema_mariae.cyst.muscle+mastersheet$Clavinema_mariae.muscle+
  mastersheet$Clavinema_mariae.flush+
  mastersheet$Clavinema_mariae.bodycav+mastersheet$Clavinema_mariae

cope.chalimis<-mastersheet$Cope_Chalimis.skin+mastersheet$Cope_Chalimis.analfin+
  mastersheet$Cope_Chalimis.caudalfin

cope.maca<-mastersheet$Cope_Maca.pectoral

copesp2<-mastersheet$Cope_unkn2.gill

acanth.dojirii<-mastersheet$Acanthochondria_dojirii.gill

ocean.pallida<-mastersheet$Oceanobdella_pallida.gill

metacercariaOS<-mastersheet$Cyst.OS.pectoral+mastersheet$Cyst.OS.muscle+
  mastersheet$Cyst.OS.pelvicfin+mastersheet$Cyst.OS.analfin+
  mastersheet$Cyst.OS.caudalfin+mastersheet$Cyst.OS.dorsal
metacercariaOS2<-decostand(x=metacercariaOS, method="pa")
as.numeric(metacercariaOS2)

trypano<-mastersheet$Trypanorhynch_larvae_cyst.intestine+
  mastersheet$Trypanorhynch_larvae_cyst.stomach

acanth.corn.echino<-mastersheet$Acanth_Corn.flush+mastersheet$Acanth_Corn.stomach+
  mastersheet$Acanth_Corn.intestine+
  mastersheet$Echinorhynchus_sp._encysted.liver+mastersheet$Echinorhynchus_lageniformis+
  mastersheet$Echinorhynchus_lageniformis.intestine+
  mastersheet$Encysted_acanthella.bodycavity+mastersheet$Encysted_acanthella.kidney+
  mastersheet$Acanth_unkn.intestine+mastersheet$Acanth_cystacanth.stomach

d.varicus<-mastersheet$Derogenes_varicus.stomach

t.lindbergi<-mastersheet$Tubulovesicula_lindbergi.stomach+
  mastersheet$Tubulovesicula_lindbergi.intestine+mastersheet$Tubulovesicula_lindbergi.gill

l.calli<-mastersheet$Lepidapedon_calli.stomach+mastersheet$Lepidapedon_calli

c.annulatus<-mastersheet$Cucullanus_annulatus.heart+
  mastersheet$Cucullanus_annulatus.intestine+mastersheet$Cucullanus_annulatus.stomach+
  mastersheet$Cucullanus_annulatus.flush+mastersheet$Cucullanus_annulatus.bodycav

contracaecum<-mastersheet$Contracaecum_larvae+
  mastersheet$Contracaecum_larvae.flush+mastersheet$Contracaecum_larvae.stomach+
  mastersheet$Contracaecum_larvae.intestine+mastersheet$Contracaecum_larvae.bodycav+
  mastersheet$Contracaecum_larvae.liver+mastersheet$Contracaecum_larvae.kidney
contracaecumsp2<-mastersheet$Contracaecum_sp2

tremsp1<-mastersheet$Trematode_sp1.stomach+mastersheet$Trematode_sp1.intestine+
  mastersheet$Trematode_sp1.flush+
  mastersheet$Trematode_sp1.gill+mastersheet$Trematode_sp1

tremsp2<-mastersheet$Trematode_sp2+mastersheet$Trematode_sp2.intestine+
  mastersheet$Trematode_sp2.stomach+
  mastersheet$Trematode_sp2.flush+mastersheet$Trematode_sp2.gill

tremsp3<-mastersheet$Trematode_sp3+mastersheet$Trematode_sp3.gill+
  mastersheet$Trematode_sp3.stomach+mastersheet$Trematode_sp3.flush

c.parophrysi<-mastersheet$Capillaria_parophrysi.stomach+
  mastersheet$Capillaria_parophrysi.intestine+mastersheet$Capillaria_parophrysi.flush

spiruridsp1<-mastersheet$Spirurid_sp1+mastersheet$Spirurid_sp1.Gill+
  mastersheet$Spirurid_sp1.Stomach+mastersheet$Spirurid_sp1.Intestine

spiruridsp2<-mastersheet$Spirurid_sp2.Stomach

oxyuridsp1<-mastersheet$Oxyurid_sp1

metacercariaA<-mastersheet$Metacercaria.gill

immaturetrem<-mastersheet$Immature_Trem.stomach

unknown.nem<-mastersheet$Unknown_nem

unknown.trem<-mastersheet$Unknown_Trem

cyst.nem<-mastersheet$Cyst.nem

hostid<-mastersheet$Fish.ID
hostsp<-mastersheet$Family.Species
long<-mastersheet$long
lat<-mastersheet$lat
year<-mastersheet$Year.Collected
tl<-mastersheet$Total.Length..mm.
sl<-mastersheet$Standard.Length..mm.
temp<-mastersheet$RaceRocks_ST_MonthYear.C

datasheet1<-data.frame(hostid, hostsp, long, lat, year, tl, sl, temp,
                       clav.mariae,cope.chalimis,cope.maca,copesp2,
                       acanth.dojirii,ocean.pallida,metacercariaOS2,metacercariaA,trypano,
                       acanth.corn.echino,d.varicus,t.lindbergi,l.calli,
                       c.annulatus,contracaecum,contracaecumsp2,tremsp1,tremsp2,tremsp3,c.parophrysi,
                       spiruridsp1,spiruridsp2,oxyuridsp1,immaturetrem, unknown.nem,unknown.trem, cyst.nem)
####For grouped abundance analyses
nematodeabun<-clav.mariae+c.annulatus+contracaecum+contracaecumsp2+
  spiruridsp1+spiruridsp2+oxyuridsp1+c.parophrysi+unknown.nem+cyst.nem

cestodeabun<-trypano

trematodeabun<-tremsp1+tremsp2+tremsp3+immaturetrem+unknown.trem+
  d.varicus+t.lindbergi+l.calli+metacercariaA

#Metacercaria kept separate
#MetacercariaOS is presence/absence
#MetacercariaA is abundance
metacercariaA
metacercariaOS2

acanthabun<-acanth.corn.echino

leechabun<-ocean.pallida

copepodabun<-cope.chalimis+cope.maca+copesp2+acanth.dojirii
datasheetGroup<-data.frame(datasheet1$hostid, datasheet1$hostsp, datasheet1$long,
                           datasheet1$lat, datasheet1$year, datasheet1$tl, datasheet1$sl,
                           datasheet1$temp, nematodeabun, cestodeabun, trematodeabun,
                           metacercariaA, metacercariaOS2,acanthabun, leechabun, copepodabun
)
)


#Jitter coordinates to remove repeats - required for autocorrelation analyses
latjitt<-jitter(datasheet1$lat, factor=0.1, amount=NULL)
longjitt<-jitter(datasheet1$long, factor=0.1, amount=0)
latjitt<-jitter(datasheetGroup$lat, factor=0.1, amount=NULL)
longjitt<-jitter(datasheetGroup$long, factor=0.1, amount=0)

datasheet<-cbind(datasheet1, longjitt, latjitt) 

colnames(datasheet)
#This is where the analysis begins!
#The below list includes all indiv parasites for the indiv models 
EngSoleParasiteList1=c( "clav.mariae","cope.chalimis","cope.maca","copesp2",
                        "acanth.dojirii","ocean.pallida","metacercariaOS2","metacercariaA","trypano",
                        "acanth.corn.echino","d.varicus","t.lindbergi","l.calli","c.annulatus",
                        "contracaecum","contracaecumsp2","tremsp1","tremsp2","tremsp3","c.parophrysi",
                        "spiruridsp1","spiruridsp2","oxyuridsp1","immaturetrem", "unknown.nem",
                        "unknown.trem","cyst.nem")
#Not every parasite we encountered occurred often enough to measure change over time. We set an arbitrary threshold oof 5% prevalence, and thus need to know which parasites are within this threshold. We'll put this info into a dataframe and save the results as its own file for later reference.
EngSolePrev=data.frame(Parasite=EngSoleParasiteList1,total=NA)
for(x in 1:length(EngSoleParasiteList1)){
  eval(parse(text=paste0('EngSolePrev$total[x]=length(datasheet$',EngSoleParasiteList1[x],'[which(datasheet$',EngSoleParasiteList1[x],'>0)])/length(datasheet$',EngSoleParasiteList1[x],')')))  
}
write.csv(EngSolePrev,'EngSolePrev_results.csv',row.names = F)
# From these results we see that 12/24 of the identified parasite species have greater than 5% prevalence. These 12 will be the ones we analyse in our individual models. So, we will now make a new list to call from. See below.

EngSoleParasiteList=c("clav.mariae","cope.chalimis","ocean.pallida", "metacercariaOS2", "metacercariaA","c.annulatus","contracaecum","tremsp1", "tremsp3","c.parophrysi","spiruridsp1")


#The group list will be used for analyses where we include all parasites by group, so that low prevalence and unknown spp. are included.
EngSoleGroupList = c("nematodeabun", "cestodeabun", "trematodeabun","acanthabun", "leechabun", "copepodabun")

#We now need to test for spatial (and temporal autocorrelation next) for each parasite. Following the entire set of spatial analyses, we need to perform an FDR adjustment of the original p-values to account for false discovery.
spatialtestsummary.parasite<-list()
spatialtestsummary.observed<-list()
spatialtestsummary.expected<-list()
spatialtestsummary.sd<-list()
spatialtestsummary.p<-list()
spatialtestsummary.fdrp<-list()
for (x in 1:length(EngSoleParasiteList)){
  parasite.data<-
    datasheet %>% 
    dplyr::select(EngSoleParasiteList[x], year, longjitt, latjitt)
  eval(parse(text=paste0('spatial.parasite<-glm(',EngSoleParasiteList[x],'~ year,
                        family="poisson", data= parasite.data)')))
  simspatial.parasite<-simulateResiduals(fittedModel = spatial.parasite)
  spatialtest<-testSpatialAutocorrelation(simulationOutput = simspatial.parasite,  x = parasite.data$longitt, y = parasite.data$latjitt)
  spatialtestsummary.parasite[[x]]<-EngSoleParasiteList[x]
  spatialtestsummary.observed[[x]]<-spatialtest$statistic[[1]]
  spatialtestsummary.expected[[x]]<-spatialtest$statistic[[2]]
  spatialtestsummary.sd[[x]]<-spatialtest$statistic[[3]]
  spatialtestsummary.p[[x]]<-spatialtest$p.value
  print(x)
}
spatialtestsummary.observed<-do.call(rbind, spatialtestsummary.observed)
spatialtestsummary.expected<-do.call(rbind, spatialtestsummary.expected)
spatialtestsummary.sd<-do.call(rbind, spatialtestsummary.sd)
spatialtestsummary.p<-do.call(rbind, spatialtestsummary.p)
spatialtestsummary.fdrp<-list(rbind, spatialtestsummary.fdrp)
spatialtestsummarystats<-do.call(cbind, list(unlist(spatialtestsummary.parasite),spatialtestsummary.observed,spatialtestsummary.expected, spatialtestsummary.sd, spatialtestsummary.p, spatialtestsummary.fdrp))
spatialtestsummarystats
spatialtestsummarystats<-as.data.frame(spatialtestsummarystats)
spatialtestsummarystats$V6=p.adjust(spatialtestsummarystats$V5,method = "fdr")
spatialtestsummarystats
colnames(spatialtestsummarystats)<-c('parasite',"obs", "exp", "sd", "p", "fdrp")
write.csv(as.matrix(spatialtestsummarystats),'spatialtest.csv')
# There was no spatial autocorr in the models before and after the fdr adjustment; going to double check metacercaria.os2 bc its binary
spatial.metos2.2<-glm(metacercariaOS2~ year, family="binomial", data= parasite.data)
simspatial.metos2.2<-simulateResiduals(fittedModel = spatial.metos2.2)
spatialtestmetos2.2<-testSpatialAutocorrelation(simulationOutput = simspatial.metos2.2,  x = parasite.data$longitt, y = parasite.data$latjitt)
spatialtestmetos2.2
#no spatial autocorr found. yay!
#We now need to test for temporal autocorrelation for each parasite. Following the entire set of spatial analyses, we need to perform an FDR adjustment of the original p-values to account for false discovery.
time<-datasheet$year
temporaltestsummary.parasite<-list()
temporaltestsummary.dw<-list()
temporaltestsummary.p<-list()
temporaltestsummary.fdrp<-list()
for (x in 1:length(EngSoleParasiteList)){
  parasitetemporal.data<-
    datasheet %>% 
    dplyr::select(EngSoleParasiteList[x], year)
  eval(parse(text=paste0('temporal.parasite<-glm(',EngSoleParasiteList[x],'~ year, family=poisson(link="log"), data= parasitetemporal.data)')))
  temporaltest<-dwtest(temporal.parasite, order.by = time, alternative = "two.sided", exact = FALSE, tol = 1e-10)
  temporaltest
  temporaltestsummary.parasite[[x]]<-EngSoleParasiteList[x]
  temporaltestsummary.dw[[x]]<-temporaltest$statistic[[1]]
  temporaltestsummary.p[[x]]<-temporaltest$p.value
  print(x)
}
temporaltestsummary.dw<-do.call(rbind, temporaltestsummary.dw)
temporaltestsummary.p<-do.call(rbind, temporaltestsummary.p)
temporaltestsummary.fdrp<-list(rbind,temporaltestsummary.fdrp)
temporaltestsummarystats<-do.call(cbind, list(unlist(temporaltestsummary.parasite),temporaltestsummary.dw,temporaltestsummary.p, temporaltestsummary.fdrp))
temporaltestsummarystats
temporaltestsummarystats<-as.data.frame(temporaltestsummarystats)
temporaltestsummarystats$V4=p.adjust(temporaltestsummarystats$V3,method = "fdr")
colnames(temporaltestsummarystats)<-c("parasite","dw","p", "fdrp") 
temporaltestsummarystats
write.csv(as.matrix(temporaltestsummarystats),'temporaltest.csv')

#from these results we can see that before and after fdrp, clav.mariae, and tremsp1 are temporally auto-correlated.
#need to check metacercaria.os2 bc its binary not poisson
temporal.metos2<-glm(metacercariaOS2~ year, family=binomial, data= parasitetemporal.data)
temporal.metos2test<-dwtest(temporal.metos2, order.by = time, alternative = "two.sided", exact = FALSE, tol = 1e-10)
temporal.metos2test
#These results are nearly the same as the above; no temp autcorr for metOS2
#Okay, so it's time to actually run the models, and I will first work on the 3 that need to account for spatial autocorr.
datasheet$sites1=paste(datasheet$lat,datasheet$long)
datasheet$logsl= log(datasheet$sl)
#clavmariae (1st to be autocorrected)
clav.mariae.mod<-glmmPQL(clav.mariae ~ I(year-1930)+ offset(logsl), random = ~1|sites1, data=datasheet, family = nbinom2(link="log"),correlation = corCAR1(), verbose = FALSE)
summary(clav.mariae.mod)  #temp corrected and not sig


clav.mariae.mod<-glmmTMB(clav.mariae~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = nbinom2(link="log")) 
summary(clav.mariae.mod) 
temporaltestclav<-dwtest(clav.mariae.mod, order.by = year, alternative = "two.sided", exact = FALSE, tol = 1e-10)


# just ran out of curiosity; not  tempcorrected and significant; this is the model for the plot!
#make the plot
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(clav.mariae.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(clav.mariae.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
clav.mariae.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Clavinema mariae") + theme(axis.title.y=element_text(face="italic")) +
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = clav.mariae/sl),
             color="#7bccc4") + #7bccc4"
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
clav.mariae.plot
pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")



#tremsp1 (rd and last to be autocorrected)
tremsp1.mod<-glmmPQL(tremsp1 ~ I(year-1930)+ offset(logsl), random = ~1|sites1, data=datasheet, family = nbinom2(link="log"),correlation = corCAR1(), verbose = FALSE)
summary(tremsp1.mod)  
tremsp1.mod2<-glmmTMB(tremsp1~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = nbinom2(link="log")) 
summary(tremsp1.mod2)
temporaltesttrem1<-dwtest(tremsp1.mod2, order.by = year, alternative = "two.sided", exact = FALSE, tol = 1e-10)

#makeplot

year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(tremsp1.mod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(tremsp1.mod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
tremsp1.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Trematoda sp.1") + 
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = tremsp1/sl),
             color="#08589e") + 
  
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
tremsp1.plot
#cope.chalimis
cope.chalimis.mod<-glmmTMB(cope.chalimis~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = nbinom2(link="log")) 
summary(cope.chalimis.mod)
rescope.cope.chamlimis = simulateResiduals(cope.chalimis.mod)
plot(rescope.cope.chamlimis)
testDispersion(rescope.cope.chamlimis)
#makeplot

year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(cope.chalimis.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(cope.chalimis.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
cope.chalimis.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Copepoda sp.") + 
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = cope.chalimis/sl),
             color="#7bccc4") + 
  
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
cope.chalimis.plot

#oceano. pallida
ocean.pallida.mod<-glmmTMB(ocean.pallida~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = nbinom2(link="log")) 
summary(ocean.pallida.mod)
rescope_ocean.pallida.mod = simulateResiduals(ocean.pallida.mod)
plot(rescope_ocean.pallida.mod)
testDispersion(rescope_ocean.pallida.mod)
#make plot
obs_pred <- predict(ocean.pallida.mod, type = "response", se.fit = TRUE)
obs_pred$year.collected <- datasheet$year
obs_pred$std_abund <- obs_pred$fit / mastersheet$sl
obs_pred <- as.data.frame(obs_pred)
obs_pred
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(ocean.pallida.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(ocean.pallida.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
ocean.pallida.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Oceanobdella pallida") + theme(axis.title.y=element_text(face="italic")) +
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = ocean.pallida/sl),
             color="#7bccc4") + 
  
  geom_line(aes(y = fit),
            data = year_pred_df) + 
geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
ocean.pallida.plot

#metacercariaOS2
metacercaria.OS2.mod2<-glmmTMB(metacercariaOS2~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = binomial) 
summary(metacercaria.OS2.mod2)
#makeplot
obs_pred <- predict(metacercaria.OS2.mod2, type = "response", se.fit = TRUE)
obs_pred$year.collected<-datasheet$year
obs_pred$std_abund <- obs_pred$fit / datasheet$sl
obs_pred<-as.data.frame(obs_pred)
obs_pred
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(metacercaria.OS2.mod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(metacercaria.OS2.mod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
metacercariaOS2.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Metacercaria sp. 2") +
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = metacercariaOS2/sl),
             color="#08589e") + 
  
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
metacercariaOS2.plot

#metacercaria.a
metacercaria.a.mod<-glmmTMB(metacercariaA ~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = nbinom2(link="log")) 
summary(metacercaria.a.mod)
rescope_metacercaria.a = simulateResiduals(metacercaria.a.mod)
plot(rescope_metacercaria.a)
testDispersion(rescope_metacercaria.a)
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(metacercaria.a.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(metacercaria.a.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
metacercaria.a.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Metacercaria sp. 1")  +
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = metacercariaA/sl),
             color="#7bccc4") + 
  
  geom_line(aes(y = fit),
            data = year_pred_df) + 
geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
metacercaria.a.plot

#c.annulatus
c.annulatus.mod<-glmmTMB(c.annulatus~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = nbinom2(link="log")) 
summary(c.annulatus.mod)
rescope.c.annulatus = simulateResiduals(c.annulatus.mod)
plot(rescope.c.annulatus)
testDispersion(rescope.c.annulatus)
obs_pred <- predict(c.annulatus.mod, type = "response", se.fit = TRUE)
df <- data.frame(matrix(unlist(obs_pred), ncol=length(obs_pred))); colnames(df) = c("fit", "se.fit")
obs_pred$year<- datahsheet$year
obs_pred$std_abund <- obs_pred$fit / datasheet$sl
obs_pred
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(c.annulatus.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(c.annulatus.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
c.annulatus.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Cucullanus annulatus") +  theme(axis.title.y=element_text(face="italic")) +
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = c.annulatus/sl),
             color="#7bccc4") + 
  
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
c.annulatus.plot

#contracaecum
contracaecum.mod<-glmmTMB(contracaecum~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = nbinom2(link="log")) 
summary(contracaecum.mod)
rescope.contracaecum.mod = simulateResiduals(contracaecum.mod)
plot(rescope.contracaecum.mod )
testDispersion(rescope.contracaecum.mod )
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(contracaecum.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(contracaecum.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
contracaecum.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Contracaecum sp.") +  #theme(axis.title.y=element_text(face="italic")) +
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = contracaecum/sl),
             color="#7bccc4") + 
  
  geom_line(aes(y = fit),
            data = year_pred_df) + 
geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
contracaecum.plot
#trem sp 3.
trem.sp3.mod<-glmmTMB(tremsp3~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = nbinom2(link="log")) 
summary(trem.sp3.mod)
rescope8 = simulateResiduals(trem.sp3.mod)
plot(rescope8)
testDispersion(rescope8)
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(trem.sp3.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(trem.sp3.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
trem.sp3.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Trematoda sp.3") +
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = tremsp3/sl),
             color="#7bccc4") + 
  
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
trem.sp3.plot

#c.parophysi
c.parophrysi.mod<-glmmTMB(c.parophrysi~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = nbinom2(link="log")) 
summary(c.parophrysi.mod)
rescope9 = simulateResiduals(c.parophrysi.mod)
plot(rescope9)
testDispersion(rescope9)
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(c.parophrysi.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(c.parophrysi.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
c.parophrysi.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Capillaria parophrysi") +  theme(axis.title.y=element_text(face="italic")) +
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = c.parophrysi/sl),
             color="#7bccc4") + 
  
  geom_line(aes(y = fit),
            data = year_pred_df) + 
geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
c.parophrysi.plot

#spiruid sp. 1
spirurid.sp1.mod<-glmmTMB(spiruridsp1~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheet, family = nbinom2(link="log")) 
summary(spirurid.sp1.mod)
rescope10 = simulateResiduals(spirurid.sp1.mod)
plot(rescope10)
testDispersion(rescope10)
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(spirurid.sp1.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(spirurid.sp1.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
spirurid.sp1.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Spirurid sp.1")  +
  xlab("Year collected") + 
  geom_point(data = datasheet,
             aes(y = spiruridsp1/sl),
             color="#7bccc4") + 
  
  geom_line(aes(y = fit),
            data = year_pred_df) + 
geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
spirurid.sp1.plot

#for the above models, I need to account for false disovery rate. see below
#now i need to adjust these p values with FDR
fdrcorrection<-read.csv("3.20.results.csv", header=TRUE)
fdrcorrection$year.fdr.p<-p.adjust(fdrcorrection$year.p,method="fdr")
write.csv(fdrcorrection,'3.20.fdr.results.csv',row.names = F)
library(cowplot)
library (patchwork)

parasiteindivplot<-plot_grid(cope.chalimis.plot,ocean.pallida.plot, tremsp1.plot, trem.sp3.plot,
                             metacercaria.a.plot,metacercariaOS2.plot, clav.mariae.plot,c.annulatus.plot,contracaecum.plot,c.parophrysi.plot, spirurid.sp1.plot)
parasiteindivplot
ggsave("parasiteindivplot.pdf")
#Now that we have looked at the inidivdial species, we will conduct the meta-analysis to understand patterns of changes with broader resolution and across transmission strategey (simple versus complex).

###Metanalysis----
library(metafor)
metadat1= data.frame(read.table('3.20.fdr.results.csv', header = T, sep=','))
metadat2=escalc(measure = "GEN", yi=year.est, sei = year.se, data = metadat1)
##ADD IN OTHER MODS HERE with interaction
metamod.null=rma.mv(yi=yi, V=vi, mods=~1,data=metadat2) #overall pattern for all parasites
metamod.group=rma.mv(yi=yi, V=vi, mods=~-1+group,data=metadat2) # patterns within parasite group\
metamod.transmission=rma.mv(yi=yi, V=vi, mods=~-1+transmission, data=metadat2)
summary(metamod.null)
summary(metamod.group)
summary(metamod.transmission)
#### plot meta-analysis 

predict.meta=data.frame(predict(metamod.group))
predict.meta$group=metadat2$group[!is.na(metadat2$year.est)]
predict.meta2=predict.meta#[c(1,2,3,5,8),]
predict.fulleffect=data.frame(predict(metamod.null))
predict.fulleffect$group='All'
predict.transmission=data.frame(predict(metamod.transmission))
predict.transmission$group=metadat2$transmission
predict.meta4=rbind(predict.meta2, predict.transmission, predict.fulleffect) #was rbind.fill
print(predict.meta4)

#Plot the predicted effects from the meta-analytic model according to parasite group and stage

metaplot1<-ggplot(aes(x=group,y=pred),data=predict.meta4)+
  geom_crossbar(aes(ymax=ci.ub,ymin=ci.lb,fill=group),width=.5,position=position_dodge(width = .5))+
  geom_hline(yintercept = 0)+ylab('Effect Size')+
  xlab('')+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(legend.position = 'none')+#c(.8,.75) would be what to replace 'none' with if a legend is desired
  scale_fill_manual(name='parasite group',limits=c("All","Direct", "Complex", "Hirudinea", "Copepoda", "Trematoda", "Nematoda"), values=c("#08589e","#2b8cbe", "#4eb3d3","#7bccc4","#a8ddb5","#ccebc5",
  "#dadaeb"), labels=c("All","Direct", "Complex", "Hirudinea", "Copepoda", "Trematoda", "Nematoda")) +
  scale_color_manual(values = c('black','black','black'))+
  theme(text=element_text(family='sans',size=20,color = 'black'),axis.text = element_text(size=20,color='black'),legend.text = element_text(size=20))+
  scale_x_discrete(limits=c("Nematoda", "Trematoda", "Copepoda","Hirudinea", "Complex", "Direct", "All"))+ ylim(c(-0.08,0.08)) +
   geom_vline(data=predict.meta4, aes(xintercept=6.6),linetype="solid", size=0.5) +
   geom_vline(data=predict.meta4, aes(xintercept=4.6),linetype="solid", size=0.5) + 
  coord_flip()
metaplot1
ggsave("metaplot1.pdf")

#The group list will be used for analyses where we include all parasites by group, so that low prevalence and unknown spp. are not a problem.
EngSoleGroupList = c("nematodeabun", "cestodeabun", "trematodeabun","acanthabun", "leechabun", "copepodabun")

#We now need to test for spatial (and temporal autocorrelation next) for each group. Following the entire set of spatial analyses, we need to perform an FDR adjustment of the original p-values to account for false discovery.
spatial.nematodeadun<-glm(nematodeabun~year, family = "poisson", data=datasheetGroup)
simspatial.nematodeabun<-simulateResiduals(fittedModel = spatial.nematodeadun)
spatialtest.nematodeabun<-testSpatialAutocorrelation(simulationOutput = simspatial.nematodeabun,  x = datasheetGroup$longitt, y = datasheetGroup$latjitt)


spatial.cestodeadun<-glm(cestodeabun~year, family = "poisson", data=datasheetGroup)
simspatial.cestodeabun<-simulateResiduals(fittedModel = spatial.cestodeadun)
spatialtest.cestodeabun<-testSpatialAutocorrelation(simulationOutput = simspatial.cestodeabun,  x = datasheetGroup$longitt, y = datasheetGroup$latjitt)

spatial.trematodeadun<-glm(trematodeabun~year, family = "poisson", data=datasheetGroup)
simspatial.trematodeabun<-simulateResiduals(fittedModel = spatial.trematodeadun)
spatialtest.trematodeabun<-testSpatialAutocorrelation(simulationOutput = simspatial.trematodeabun,  x = datasheetGroup$longitt, y = datasheetGroup$latjitt)


spatial.acanthadun<-glm(acanthabun~year, family = "poisson", data=datasheetGroup)
simspatial.acanthabun<-simulateResiduals(fittedModel = spatial.acanthadun)
spatialtest.acanthabun<-testSpatialAutocorrelation(simulationOutput = simspatial.acanthabun,  x = datasheetGroup$longitt, y = datasheetGroup$latjitt)

spatial.leechadun<-glm(leechabun~year, family = "poisson", data=datasheetGroup)
simspatial.leechabun<-simulateResiduals(fittedModel = spatial.leechadun)
spatialtest.leechabun<-testSpatialAutocorrelation(simulationOutput = simspatial.leechabun,  x = datasheetGroup$longitt, y = datasheetGroup$latjitt)

spatial.copepodadun<-glm(copepodabun~year, family = "poisson", data=datasheetGroup)
simspatial.copepodabun<-simulateResiduals(fittedModel = spatial.copepodadun)
spatialtest.copepodabun<-testSpatialAutocorrelation(simulationOutput = simspatial.copepodabun,  x = datasheetGroup$longitt, y = datasheetGroup$latjitt)
#none of the above spatial tests were significant, so runnings FDRs is a moot point.
#lets test the group models for temporal autocorr

timegroup<-datasheetGroup$year
temporal.nematodeabun<-glm(nematodeabun~ year, family=poisson(link="log"), data= datasheetGroup)
temporaltest.nematodeabun<-dwtest(temporal.nematodeabun, order.by = timegroup, alternative = "two.sided", exact = FALSE, tol = 1e-10) #temporal autocorrelation present


temporal.cestodeabun<-glm(cestodeabun~ year, family=poisson(link="log"), data= datasheetGroup)
temporaltest.cestodeabun<-dwtest(temporal.cestodeabun, order.by = timegroup, alternative = "two.sided", exact = FALSE, tol = 1e-10) #temporal autocorrelation not present


temporal.trematodeabun<-glm(trematodeabun~ year, family=poisson(link="log"), data= datasheetGroup)
temporaltest.trematodeabun<-dwtest(temporal.trematodeabun, order.by = timegroup, alternative = "two.sided", exact = FALSE, tol = 1e-10)  #temp autcorr not present


temporal.acanthabun<-glm(acanthabun~ year, family=poisson(link="log"), data= datasheetGroup)
temporaltest.acanthabun<-dwtest(temporal.acanthabun, order.by = timegroup, alternative = "two.sided", exact = FALSE, tol = 1e-10) #temp autocorr not present

temporal.leechabun<-glm(leechabun~ year, family=poisson(link="log"), data= datasheetGroup)
temporaltest.leechabun<-dwtest(temporal.leechabun, order.by = timegroup, alternative = "two.sided", exact = FALSE, tol = 1e-10) # temp autrcorr not present

temporal.copepodabun<-glm(copepodabun~ year, family=poisson(link="log"), data= datasheetGroup)
temporaltest.copepodabun<-dwtest(temporal.copepodabun, order.by = timegroup, alternative = "two.sided", exact = FALSE, tol = 1e-10)# temp autrcorr not present
#nematode abundace is temporally autocorrelated, so I will fdr-adjust; likely still sig.
grouppvalues<-c(1.31E-08,
                0.8704,
                0.3666,
                0.1623,
                0.6192,
                0.7568)
grouppvaluesfdr<-p.adjust(grouppvalues, method ="BH", n=length(grouppvalues))
#nematode abundance is still temporally autocorrelated, so I will account for that in the model.

datasheetGroup$sites1=paste(datasheet$lat,datasheet$long)
datasheetGroup$logsl= log(datasheet$sl)
nematode.abun.mod<-glmmPQL(nematodeabun ~ I(year-1930)+ offset(logsl), random = ~1|sites1, data=datasheetGroup, family = nbinom2(link="log"),correlation = corCAR1(), verbose = FALSE)
summary(nematode.abun.mod)  #temp corrected and not sig
nematode.abun.mod2<-glmmTMB(nematodeabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = nbinom2(link="log")) 
summary(nematode.abun.mod2) # just ran out of curiosity; not  tempcorrected and significant
#make the plot
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(nematode.abun.mod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(nematode.abun.mod2, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
nematodeabun.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Nematoda") + #theme(axis.title.y=element_text(face="italic")) +
  xlab("Year collected") + 
  geom_point(data = datasheetGroup,
             aes(y = nematodeabun/sl),
             color="#7bccc4") + 
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))
nematodeabun.plot

#cestodes
cestode.abun.mod<-glmmTMB(cestodeabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = nbinom2(link="log")) #issues with convergence
summary(cestode.abun.mod)
cestodedispersion= simulateResiduals(cestode.abun.mod)
plot(cestodedispersion)
testDispersion(cestodedispersion) # overdispersed!
cestode.abun.mod.pois<-glmmTMB(cestodeabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = "poisson") #the poisson model AIC is lower by ~2, the results are identical; reporting this
summary(cestode.abun.mod.pois)
#make the plot
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(cestode.abun.mod.pois, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(cestode.abun.mod.pois, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
cestodeabun.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Cestoda")  + ylim(0,1) +#added this in to fix odd axis scale but needs to be addressed further 
  #xlab("Year collected") 
  geom_point(data = datasheetGroup,
             aes(y = cestodeabun/sl),
             color="#7bccc4") + 
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black')) +xlab("")
cestodeabun.plot

#trematode abund
trematode.abun.mod<-glmmTMB(trematodeabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = nbinom2(link="log")) 
summary(trematode.abun.mod)
trematodedispersion= simulateResiduals(trematode.abun.mod)
plot(trematodedispersion)
testDispersion(trematodedispersion) # NOT Overdispsed
#trematode.abun.mod.pois<-glmmTMB(trematodeabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = "poisson") #the poisson model AIC is lower by ~350
#summary(trematode.abun.mod.pois) #lower AIC
#make the plot
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(trematode.abun.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(trematode.abun.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
trematodeabun.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Trematoda")  +
  #xlab("Year collected") + 
  geom_point(data = datasheetGroup,
             aes(y = trematodeabun/sl),
             color="#7bccc4") + 
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black')) + xlab("")
trematodeabun.plot

#acanthabun
acanth.abun.mod<-glmmTMB(acanthabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = nbinom2(link="log")) 
summary(acanth.abun.mod)
acanthabundispersion= simulateResiduals(acanth.abun.mod)
plot(acanthabundispersion)
testDispersion(acanthabundispersion) # NOT Overdispsed
#acanthabun.abun.mod.pois<-glmmTMB(acanthabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = "poisson") #AIC 3 higher.
#summary(acanthabun.abun.mod.pois)
#make the plot
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(acanth.abun.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(acanth.abun.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
acanthabun.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Acanthocephala")  +
  #xlab("Year collected") + 
  geom_point(data = datasheetGroup,
             aes(y = acanthabun/sl),
             color="#08589e") + 
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black'))+ xlab("")
acanthabun.plot

#leech abund
leech.abun.mod<-glmmTMB(leechabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = nbinom2(link="log")) 
summary(leech.abun.mod)
leechabundispersion= simulateResiduals(leech.abun.mod)
plot(leechabundispersion)
testDispersion(leechabundispersion) # NOT Overdispsed
#leechabun.abun.mod.pois<-glmmTMB(leechabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = "poisson") #AIC 6 higher.
#summary(leechabun.abun.mod.pois)
#make the plot
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(leech.abun.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(leech.abun.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
leechabun.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Hiruduina")  +
  #xlab("Year collected") + 
  geom_point(data = datasheetGroup,
             aes(y = leechabun/sl),
             color="#7bccc4") + 
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black')) + xlab("")
leechabun.plot

#copepodabun
copepod.abun.mod<-glmmTMB(copepodabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = nbinom2(link="log")) 
summary(copepod.abun.mod)
copepodabundispersion= simulateResiduals(copepod.abun.mod)
plot(copepodabundispersion)
testDispersion(copepodabundispersion) # NOT Overdispsed
#copepodabun.abun.mod.pois<-glmmTMB(copepodabun~ I(year-1930)+ offset(logsl) + (1|sites1),data=datasheetGroup, family = "poisson") #AIC 12 higher.
#summary(copepodabun.abun.mod.pois)
#make the plot
year_pred_df <- data.frame(year = 1930:2019,
                           logsl = 0, # log(1) = 0,
                           sites1 = 100)
year_pred <- predict(copepod.abun.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "response")
year_pred_df <- cbind(year_pred_df, year_pred)
year_quantiles <- predict(copepod.abun.mod, year_pred_df, allow.new.levels = TRUE, se.fit = TRUE, type = "link")
year_quantiles$upper_log <- year_quantiles$fit + 1.96 * year_quantiles$se.fit
year_quantiles$lower_log <- year_quantiles$fit - 1.96 * year_quantiles$se.fit
year_quantiles$upper <- exp(year_quantiles$upper_log)
year_quantiles$lower <- exp(year_quantiles$lower_log)
year_quantiles <- as.data.frame(year_quantiles)
year_quantiles$year= 1930:2019
copepodabun.plot<-ggplot(mapping = aes(x = year)) +
  theme_classic() +
  ylab("Copepoda")  +
  #xlab("Year collected")  
  geom_point(data = datasheetGroup,
             aes(y = copepodabun/sl),
             color="#7bccc4") + 
  geom_line(aes(y = fit),
            data = year_pred_df) + 
  geom_ribbon(data = year_quantiles, aes(ymin = lower, ymax = upper ), fill="gray", alpha=0.4 ) + theme(text=element_text(size=12,color='black')) + xlab("")
copepodabun.plot
#the result of the groupabun models now need to be fdr-adjusted
groupmodels<-c(0.0335,
               0.60772,
               0.116,
               0.000436,
               0.159,
               0.354)
groupmodelsfdr<-p.adjust(groupmodels, method ="BH", n=length(groupmodels))
library(cowplot)
parasitegroupplot<-plot_grid(copepodabun.plot,acanthabun.plot, leechabun.plot,trematodeabun.plot,nematodeabun.plot,cestodeabun.plot)
parasitegroupplot
ggsave("parasitegroupplot.pdf")
#lets work on some descriptive stats
#What is the total # of indiv parasites found? What percent does each group make up?
totalsum<-sum(datasheetGroup$nematodeabun, datasheetGroup$cestodeabun, datasheetGroup$trematodeabun, datasheetGroup$acanthabun, datasheetGroup$leechabun, datasheetGroup$copepodabun)
#nematode proportion
sum(datasheetGroup$nematodeabun)/totalsum
#cestode proportion= 
sum(datasheetGroup$cestodeabun)/totalsum
#trematode proportion = 
sum(datasheetGroup$trematodeabun)/totalsum
#acanth proportion = 
sum(datasheetGroup$acanthabun)/totalsum
#leech prop = 
sum(datasheetGroup$leechabun)/totalsum
#copepod prop = 
sum(datasheetGroup$copepodabun)/totalsum
write.csv(datasheetGroup, "datasheetGroup.csv")
mean(datasheet$tl)
sd(datasheet$tl)
