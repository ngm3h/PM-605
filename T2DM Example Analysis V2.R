# ===== Clear R enviroment ========================================================================
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()

lapply(names(sessionInfo()$otherPkgs), function(pkgs)
  detach(
    paste0('package:', pkgs),
    character.only = T,
    unload = T,
    force = T
  ))

.rs.restartR()


# ===== Import Libraries ==========================================================================

library(dplyr)
library(ggplot2)
library(meta) # for metagen, update.meta, forest, and funnel

# ===== Sent Save Loc =============================================================================

save.loc = "E:/PM-605/Week5-Feb.22/"
data.loc = "E:/PM-605/Week5-Feb.22/"


# ===== In-Class Examples to Show Plots ===========================================================

# ===== * Load Data =====
dmbc_flw_rm_long_sexcomb_out = paste0(data.loc, 
                                      "PM605DiabetesAndBlad_DATA_2019-05-03_1706_CLEANED.rds") %>% 
  readRDS()


# ===== * Look at Data =====
dmbc_flw_rm_long_sexcomb_out %>% glimpse()

# first_auth: first author
# datepub: publication date
# location: study location
# ratio_est: ratio estimte (eg, HR, OR)
# ratio_est_l95: lower 95% CI of ratio estimate
# ratio_est_u95: upper 95% CI of ratio estimate
# study_design: study design
dmbc_flw_rm_long_sexcomb_out %>% group_by(study_design) %>% tally()

# datepub_plot: Publication date variable used for plotting (since some studies are 
# published in the same year, but different months, this var adds decimal to the 
# year var specifying time of year)



# ===== * Prep Data =====

ma_sex_comb.test = 
  dmbc_flw_rm_long_sexcomb_out %>% 
  mutate(
    lnes = log(ratio_est),
    lnes_l = log(ratio_est_l95),
    lnes_u = log(ratio_est_u95),
    selnes = (lnes_u-lnes_l)/3.919928,
    Study = paste0(first_auth," et al., " , datepub),
    cases = affected_exp + affected_unexp,
    
    # order design variable show it shows correctly in the plot
    Design = study_design %>% forcats::fct_relevel(c("1. Pro Cht","2. Pop CC","3. Ret Cht",
                                                           "4. SIR","5. Hosp CC"))
  )



# ===== * Run Sex-Combined Meta-Analysis =====
test.ma <- metagen(lnes, selnes, data = ma_sex_comb.test, studlab = paste(Study),
                comb.fixed = FALSE, comb.random = TRUE,  method.tau = "DL",
                hakn = F, prediction = FALSE, sm = "HR")
test.ma

# ===== * Update the meta object to stratify by design =====
test.ma.subgroup = update.meta(test.ma, 
                            byvar=Design, 
                            comb.random = TRUE, 
                            comb.fixed = FALSE)
test.ma.subgroup


# ===== * Plot meta-analysis =====

png(paste0(save.loc,"PDAC_airpol_forest_v2.png"), 
    height = 13, width = 6, units='in', res=600)

forest(test.ma.subgroup,
       
       sortvar = datepub_plot,
       
       rightcols = c("effect.ci", "w.random"),
       rightlabs = c("Ratio Est [95% CI]", "Weight"),
       
       leftcols = c("studlab"),
       leftlab = c("Author, Year"),
       
       just.addcols = "center",
       spacing = 1,
       
       print.tau2 = FALSE,
       print.I2.ci = FALSE,
       
       ref=1,
       col.square = "grey30",
       resid.hetstat=F,
       fontsize = 8)

dev.off()



# ===== Run Sex-Combined Meta-Analysis (high quality only) ========================================

ma_sex_comb.test.highqual = 
  ma_sex_comb.test %>% 
  filter(study_design %in% c("1. Pro Cht", "2. Pop CC"))


# ===== * Run Sex-Combined Meta-Analysis =====
test.ma.hq <- metagen(lnes, selnes, data = ma_sex_comb.test.highqual, studlab = paste(Study),
                      comb.fixed = FALSE, comb.random = TRUE,  method.tau = "DL",
                      hakn = F, prediction = FALSE, sm = "HR")
test.ma.hq

# ===== * Update the meta object to stratify by design =====
test.ma.subgroup.hq = update.meta(test.ma.hq,
                                  byvar=Design, 
                                  comb.random = TRUE, 
                                  comb.fixed = FALSE)
test.ma.subgroup.hq


# ===== * Plot meta-analysis =====

png(paste0(save.loc,"PDAC_airpol_forest_hq_v2.png"), 
    height = 13, width = 6, units='in', res=600)

# create plot from "meta" package, which takes meta object
forest(test.ma.subgroup.hq,
       
       sortvar = datepub_plot,
       
       rightcols = c("effect.ci", "w.random"), # attributes from metagen object
       rightlabs = c("Ratio Est [95% CI]", "Weight"), # names you want for objects on plot
       
       leftcols = c("studlab"), # Information you want on left column
       leftlab = c("Author, Year"),
       
       just.addcols = "center", # justification
       spacing = 1, # spacing option
       
       print.tau2 = FALSE, # do not print tau parameter
       print.I2.ci = FALSE, # do not print CI for I2 value
       
       xlim = c(.8, 3),
       at = seq(.8,3,.2),
       
       ref=1, # set reference line
       col.square = "grey30",
       resid.hetstat=F, # do not print residual heterogenity statistics
       fontsize = 8)

dev.off()


# ===== * Funnel Plot =====

# same meta-analysis code from earlier, but removeing the "sm = "HR" option of back-transforming
test.ma.hq.f <- metagen(lnes, selnes, data = ma_sex_comb.test.highqual, studlab = paste(Study),
                      comb.fixed = FALSE, comb.random = TRUE,  method.tau = "DL",
                      hakn = F, prediction = FALSE)
test.ma.hq.f


# funnel plot object from meta package, takes meta object
funnel(test.ma.hq.f, xlab="Ratio Estimate")



# ===== * Cumulative Plot =====

# run a random-effects meta-analysis (same as earler) but this time to get a "metafor" object
cuml_meta = metafor::rma(yi=lnes, sei=selnes, data=ma_sex_comb.test.highqual, slab = Study)

# pass meta-analysis object to "cuml" function, ordering by weight
cuml_meta_plt = 
  metafor::cumul(cuml_meta, order=order(-ma_sex_comb.test.highqual$selnes)) %>% 
  # convert results to dataframe
  as.data.frame() %>% 
  # generate ratio-estimtes
  mutate(Study = rownames(.),
         OR = exp(estimate),
         LOR = exp(estimate - 1.96*se),
         UOR = exp(estimate + 1.96*se)) %>% 
  
  # merge in wieghts
  left_join(ma_sex_comb.test.highqual[,c("Study","selnes")]) %>% 
  
  # order "Study" variable by weight
  mutate(Study = Study %>% forcats::fct_reorder(selnes))
  
# plot results

ggplot(data=cuml_meta_plt, aes(y=Study, x=OR, xmin=LOR, xmax=UOR))+ 
  
  geom_vline(xintercept=exp(test.ma.hq$TE.random), 
             color = "darkgrey", linetype="dashed", size = 1) +
  
  geom_point() + 
  geom_errorbarh(height=.1)+
  theme_light() +
  xlab("Ratio Estimate")


# ===== * Leave-one-out plot =====

# function to run the meta-analysis leaving out the study passed to the function

metainf = metainf(test.ma.hq, pooled="random", 
                  sortvar=ma_sex_comb.test.highqual$datepub_plot)

forest(metainf,
       
       rightcols = c("effect.ci"), # attributes from metagen object
       rightlabs = c("Ratio Est [95% CI]"), # names you want for objects on plot
       
       leftcols = c("studlab"), # Information you want on left column
       leftlab = c("Author, Year"),
       
       just.addcols = "center", # justification
       spacing = 1, # spacing option
       
       print.tau2 = FALSE, # do not print tau parameter
       print.I2.ci = FALSE, # do not print CI for I2 value
       
       xlim = c(1, 1.5),
       at = seq(1,1.5,.1),
       
       ref=1, # set reference line
       col.square = "grey30",
       resid.hetstat=F, # do not print residual heterogenity statistics
       fontsize = 8)

if (!require(devtools)) install.packages('devtools')
library(devtools)
install_version("meta", version = "4.13-0", repos = "http://cran.us.r-project.org")
# when prompted, skip updates
























