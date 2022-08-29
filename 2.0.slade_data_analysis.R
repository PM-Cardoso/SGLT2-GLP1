####################
## Description:
##  - This file includes a short descriptive analyses of the some of the dataset quirks.
####################


# Used in slade to ensure the library being used is my personal library
.libPaths(.libPaths()[c(2,1,3)])


## bartMachine: increase memery usage to 40gb of RAM
options(java.parameters = "-Xmx40g")


library(tidyverse)
library(bartMachine)
library(ggplot2)
library(ggthemes)
library(rms)

## path to output folder
output_path <- "Samples"
## make directory for outputs

dir.create(output_path)

output_path <- "Samples/SGLT2-GLP1"

## make directory for outputs
dir.create(output_path)

## make directory for outputs
dir.create(paste0(output_path,"/Data_Analysis"))

## make directory for outputs
dir.create("Plots")

###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# name: final.all
load(paste0(output_path, "/datasets/cprd_19_sglt2glp1_allcohort.Rda"))


###########################
###### Drugs prescribed by year
###########################

subtype_year <- final.all %>%
  select(glp1subtype, sglt2subtype, yrdrugstart) %>%
  unite("Subtype", glp1subtype:sglt2subtype, remove = TRUE, sep= "") %>%
  table() %>%
  prop.table(margin = 2) %>%
  as.data.frame() %>%
  mutate(Therapy = ifelse(Subtype == "cana"|Subtype == "dapa"|Subtype == "empa", "SGLT2", "GLP1"),
         Freq = 100*Freq) %>%
  mutate(Freq = ifelse(Freq == 0, NA, Freq))


plot_subtype_year <- subtype_year %>%
  ggplot() +
  theme_bw() +
  geom_point(aes(x = yrdrugstart, y = Freq, colour = Subtype)) +
  geom_path(aes(x = yrdrugstart, y = Freq, colour = Subtype, group = Subtype, linetype = Therapy)) +
  ylab("% total new drug prescriptions") +
  ggtitle("Drugs prescribed by year") +
  theme(legend.position = "right",
        axis.title.x = element_blank())


plot_subtype_year_facet <- subtype_year %>%
  ggplot() +
  theme_bw() +
  geom_point(aes(x = yrdrugstart, y = Freq, colour = Subtype)) +
  geom_path(aes(x = yrdrugstart, y = Freq, colour = Subtype, group = Subtype, linetype = Therapy)) +
  facet_wrap(~Therapy, ncol = 1) + 
  ylab("% total new drug prescriptions") +
  ggtitle("Drugs prescribed by year") +
  theme(legend.position = "right",
        axis.title.x = element_blank())




##########################
###### Change in HbA1c by subtype
##########################



# hba1c_change <- final.all %>%
#   select(glp1subtype, sglt2subtype, posthba1c_final, prehba1cmmol) %>%
#   unite("Subtype", glp1subtype:sglt2subtype, remove = TRUE, sep= "") %>%
#   mutate(HbA1c = posthba1c_final - prehba1cmmol) %>%
#   select(-posthba1c_final, -prehba1cmmol) %>%
#   group_by(Subtype) %>%
#   mutate(mean = mean(HbA1c),
#          `5%` = quantile(HbA1c, probs = c(0.05)),
#          `95%` = quantile(HbA1c, probs = c(0.95))) %>%
#   select(-HbA1c) %>%
#   ungroup() %>%
#   unique() %>%
#   mutate(Therapy = ifelse(Subtype == "cana"|Subtype == "dapa"|Subtype == "empa", "SGLT2", "GLP1"))
#   
# 
# plot_hba1c_change <- hba1c_change %>%
#   ggplot() +
#   theme_bw() +
#   geom_point(aes(x = mean, y = Subtype, colour = Therapy)) +
#   geom_errorbar(aes(xmin = `5%`, xmax = `95%`, y = Subtype, colour = Therapy), width = 0.2) +
#   xlab("HbA1c change (mmol/mol)") +
#   ylab("Therapy subtype") +
#   ggtitle("Observed HbA1c change by Therapy subtype") +
#   theme(legend.position = "right")



##########################
###### Change in HbA1c by subtype over time
##########################


baseline_dataset <- final.all %>%
  unite("Subtype", glp1subtype:sglt2subtype, remove = TRUE, sep= "") 



m1 <- ols(posthba1c_final~prehba1cmmol+rcs(yrdrugstart,3)*Subtype + drugline + ncurrtx + rcs(hba1cmonth,3)*Subtype + rcs(agetx,3)*Subtype,data=baseline_dataset,x=TRUE,y=TRUE)

m1.p <- Predict(m1,prehba1cmmol=median(baseline_dataset$prehba1cmmol),Subtype=c("cana", "dapa", "dulaglutide", "empa", "exenatide", "liraglutide", "lixisenatide"), yrdrugstart=c(seq(2007,2019)), agetx = mean(final.all$agetx), ncurrtx=2, drugline=c(2,3,4,5), hba1cmonth=6) %>%
  as.data.frame() %>% 
  mutate(Response = yhat-prehba1cmmol,
         ci.l = lower - prehba1cmmol,
         ci.u = upper - prehba1cmmol) %>%
  mutate(Response = ifelse(Subtype == "cana" & yrdrugstart < 2014, 
                           NA,
                           ifelse(Subtype == "dapa" & yrdrugstart < 2013,
                                  NA, 
                                  ifelse(Subtype == "empa" & yrdrugstart < 2014, 
                                         NA,
                                         ifelse(Subtype == "dulaglutide" & yrdrugstart < 2015,
                                                NA, 
                                                ifelse(Subtype == "exenatide",
                                                       Response,
                                                       ifelse(Subtype == "liraglutide" & yrdrugstart < 2009,
                                                              NA, 
                                                              ifelse(Subtype == "lixisenatide" & yrdrugstart < 2013,
                                                                     NA, 
                                                                     Response))))))),
         ci.l = ifelse(Subtype == "cana" & yrdrugstart < 2014, 
                       NA,
                       ifelse(Subtype == "dapa" & yrdrugstart < 2013,
                              NA, 
                              ifelse(Subtype == "empa" & yrdrugstart < 2014, 
                                     NA,
                                     ifelse(Subtype == "dulaglutide" & yrdrugstart < 2015,
                                            NA, 
                                            ifelse(Subtype == "exenatide",
                                                   ci.l,
                                                   ifelse(Subtype == "liraglutide" & yrdrugstart < 2009,
                                                          NA, 
                                                          ifelse(Subtype == "lixisenatide" & yrdrugstart < 2013,
                                                                 NA, 
                                                                 ci.l))))))),
         ci.u = ifelse(Subtype == "cana" & yrdrugstart < 2014, 
                       NA,
                       ifelse(Subtype == "dapa" & yrdrugstart < 2013,
                              NA, 
                              ifelse(Subtype == "empa" & yrdrugstart < 2014, 
                                     NA,
                                     ifelse(Subtype == "dulaglutide" & yrdrugstart < 2015,
                                            NA, 
                                            ifelse(Subtype == "exenatide",
                                                   ci.u,
                                                   ifelse(Subtype == "liraglutide" & yrdrugstart < 2009,
                                                          NA, 
                                                          ifelse(Subtype == "lixisenatide" & yrdrugstart < 2013,
                                                                 NA, 
                                                                 ci.u)))))))) %>%
  mutate(Therapy = ifelse(Subtype == "cana" | Subtype == "dapa" | Subtype == "empa", "SGLT2", "GLP1"))


plot_predicted_subtype <- m1.p %>%
  filter(drugline == 3) %>%
  ggplot() +
  theme_bw() +
  geom_line(aes(x = yrdrugstart, y = Response, colour = Subtype, linetype = Therapy)) +
  facet_wrap(~Therapy, ncol = 1) +
  scale_x_continuous("", breaks = seq(min(m1.p$yrdrugstart), max(m1.p$yrdrugstart), by = 1)) +
  ylab("Predicted Response (mmol/mol)") +
  ggtitle("Predicted HbA1c change, baseline adjusted")


plot_predicted_subtype_drugline <- m1.p %>%
  ggplot() +
  theme_bw() +
  geom_line(aes(x = yrdrugstart, y = Response, colour = Subtype, linetype = Therapy)) +
  facet_wrap(drugline~Therapy, ncol = 2) +
  # scale_x_continuous("", breaks = seq(min(m1.p$yrdrugstart), max(m1.p$yrdrugstart), by = 1)) +
  ylab("Predicted Response (mmol/mol)") +
  ggtitle("Predicted HbA1c change by drugline, baseline adjusted")






# hba1c_change_year <- final.all %>%
#   select(glp1subtype, sglt2subtype, posthba1c_final, prehba1cmmol, yrdrugstart) %>%
#   unite("Subtype", glp1subtype:sglt2subtype, remove = TRUE, sep= "") %>%
#   mutate(HbA1c = posthba1c_final - prehba1cmmol) %>%
#   select(-posthba1c_final, -prehba1cmmol) %>%
#   mutate(Therapy = ifelse(Subtype == "cana"|Subtype == "dapa"|Subtype == "empa", "SGLT2", "GLP1"))
# 
# 
# plot_hba1c_change_year <- hba1c_change_year %>%
#   ggplot() +
#   theme_bw() +
#   geom_smooth(aes(x = yrdrugstart, y = HbA1c, linetype = Therapy, colour = Subtype), formula = y ~ s(x, bs = "cs", k=5)) +
#   ylab("HbA1c change (mmol/mol)") +
#   ggtitle("Observed HbA1c change by Therapy subtype per year") +
#   scale_x_continuous("", breaks = seq(min(hba1c_change_year$yrdrugstart), max(hba1c_change_year$yrdrugstart), by = 1)) +
#   theme(legend.position = "right")
# 
# 
# plot_hba1c_change_year_facet <- hba1c_change_year %>%
#   ggplot() +
#   theme_bw() +
#   geom_smooth(aes(x = yrdrugstart, y = HbA1c, linetype = Therapy, colour = Subtype), formula = y ~ s(x, bs = "cs", k=5)) +
#   facet_wrap(~Therapy, ncol = 1) +
#   ylab("HbA1c change (mmol/mol)") +
#   ggtitle("Observed HbA1c change by Therapy subtype per year") +
#   scale_x_continuous("", breaks = seq(min(hba1c_change_year$yrdrugstart), max(hba1c_change_year$yrdrugstart), by = 1)) +
#   theme(legend.position = "right")


##########################
###### Initial HbA1c by subtype over time
##########################

initial_hba1c <- final.all %>%
  select(glp1subtype, sglt2subtype, prehba1cmmol, yrdrugstart) %>%
  unite("Subtype", glp1subtype:sglt2subtype, remove = TRUE, sep= "") %>%
  group_by(Subtype, yrdrugstart) %>%
  mutate(HbA1c = mean(prehba1cmmol)) %>%
  ungroup() %>%
  select(-prehba1cmmol) %>%
  mutate(Therapy = ifelse(Subtype == "cana" | Subtype == "dapa" | Subtype == "empa", "SGLT2", "GLP1")) %>%
  unique()

plot_initial_hba1c <- initial_hba1c %>%
  ggplot() +
  theme_bw() +
  geom_line(aes(x = yrdrugstart, y = HbA1c, colour = Subtype, linetype = Therapy)) +
  facet_wrap(~Therapy, ncol = 1) +
  scale_x_continuous("", breaks = seq(min(m1.p$yrdrugstart), max(m1.p$yrdrugstart), by = 1)) +
  ylab("Average initial HbA1c (mmol/mol)") +
  ggtitle("Initial HbA1c by subtype, per year")


##########################
###### Drugline, subtype, year of prescription
##########################


drugline_subtype_year_dataset <- final.all %>%
  select(glp1subtype, sglt2subtype, drugline, yrdrugstart) %>%
  unite("Subtype", glp1subtype:sglt2subtype, remove = TRUE, sep= "")


for (i in 1:length(unique(drugline_subtype_year_dataset$Subtype))) {
  if (i == 1) {
    drugline_subtype_year <- final.all %>%
      select(glp1subtype, sglt2subtype, drugline, yrdrugstart) %>%
      unite("Subtype", glp1subtype:sglt2subtype, remove = TRUE, sep= "") %>%
      filter(Subtype == unique(drugline_subtype_year_dataset$Subtype)[i]) %>%
      select(-Subtype) %>%
      table() %>%
      prop.table(margin = 2) %>%
      as.data.frame() %>%
      cbind(Subtype = unique(drugline_subtype_year_dataset$Subtype)[i]) %>%
      mutate(Therapy = ifelse(Subtype == "cana"|Subtype == "dapa"|Subtype == "empa", "SGLT2", "GLP1"),
             Freq = 100*Freq) %>%
      mutate(Freq = ifelse(Freq == 0, NA, Freq))
  } else {
    drugline_subtype_year <- rbind(
      drugline_subtype_year,
      final.all %>%
        select(glp1subtype, sglt2subtype, drugline, yrdrugstart) %>%
        unite("Subtype", glp1subtype:sglt2subtype, remove = TRUE, sep= "") %>%
        filter(Subtype == unique(drugline_subtype_year_dataset$Subtype)[i]) %>%
        select(-Subtype) %>%
        table() %>%
        prop.table(margin = 2) %>%
        as.data.frame() %>%
        cbind(Subtype = unique(drugline_subtype_year_dataset$Subtype)[i]) %>%
        mutate(Therapy = ifelse(Subtype == "cana"|Subtype == "dapa"|Subtype == "empa", "SGLT2", "GLP1"),
               Freq = 100*Freq) %>%
        mutate(Freq = ifelse(Freq == 0, NA, Freq))
    )
  }
}


plot_drugline_subtype_year <- drugline_subtype_year %>%
  ggplot() +
  theme_bw() +
  geom_point(aes(x = yrdrugstart, y = Freq, colour = drugline)) +
  geom_path(aes(x = yrdrugstart, y = Freq, colour = drugline, group = drugline, linetype = Therapy)) +
  facet_wrap(~Subtype, ncol = 1) +
  ylab("% drugs previously prescribed") +
  ggtitle("Drugs previously prescribed by Therapy subtype") +
  theme(legend.position = "right",
        axis.title.x = element_blank())


##########################
###### Association between baseline HbA1c and 6 month HbA1c response (outcome HbA1c – baseline HbA1c)
##########################


hba1c.association <- final.all %>%
  mutate(hba1c.change = posthba1c_final - prehba1cmmol) %>%
  select(prehba1cmmol, hba1c.change, drugclass)

plot.hba1c.histogram <- hba1c.association %>%
  ggplot() +
    theme_void() +
    geom_histogram(aes(x = prehba1cmmol))

plot.hba1c.association <- hba1c.association %>%
  ggplot() +
    theme_bw() +
    stat_smooth(aes(x = prehba1cmmol, y = hba1c.change, colour = drugclass), method = "auto") +
    xlab("Baseline HbA1c (mmol/mol)") + ylab("HbA1c Response (mmol/mol)") +
    ggtitle("Association between baseline HbA1c and HbA1c response") +
    scale_colour_manual(values=c("#998ec3","#f1a340")) +
    theme(legend.title = element_blank(),
          legend.position = c(0.80, 0.80))

# Add marginal histograms

plot.hba1c.association.marginals <- cowplot::plot_grid(
  
  plot.hba1c.association
  
  ,
  
  cowplot::plot_grid(
    
    ggplot() +
      theme_void()
    
    ,
    
    plot.hba1c.histogram
    
    , ncol = 2, nrow = 1, rel_widths = c(0.06, 1)
    
  )
  
  , ncol = 1, nrow = 2, rel_heights = c(0.85, 0.15)
  
)




###############################################################################
###############################################################################
###############################################################################




pdf(file = "Plots/2.0.analysis.pdf")
plot_subtype_year
plot_subtype_year_facet
plot_predicted_subtype
plot_predicted_subtype_drugline
plot_initial_hba1c
plot_drugline_subtype_year
plot.hba1c.association.marginals
dev.off()







