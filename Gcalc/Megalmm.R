setwd("~/Google Drive/Alas Drosophila/FULL DIM WD_REVISION/Gcalc/")
library(geomorph)
# library(MCMCglmm)
library(MegaLMM)
library(psych)
library(reshape2)
library(plyr)
library(cowplot)
library(dplyr)
library(readxl)
library(ggplot2)
source('~/Dropbox/__R/Rfun/ellipse.Rfun.R')

load("../GPA/GPA.Rdata")

GPA$data

morph<-data.frame(GPA$data[,1:30])

# morph.s<-scale(morph,center = F)
# data<- data.frame(morph,classifier)
# data$Sexo<-factor(data$Sexo)
# 
# P_m<-var(lm(as.matrix(morph.s)~classifier$Poblacion,subset = classifier$Sexo=="M")$res)
# E_m<-var(lm(as.matrix(morph.s)~classifier$Linea2,subset = classifier$Sexo=="M")$res)
# P_f<-var(lm(as.matrix(morph.s)~classifier$Poblacion,subset = classifier$Sexo=="H")$res)
# E_f<-var(lm(as.matrix(morph.s)~classifier$Linea2,subset = classifier$Sexo=="H")$res)
# G_m<-P_m-E_m
# G_f<-P_f-E_f

run_parameters = MegaLMM_control(
  max_NA_groups = 3,
  scale_Y = TRUE,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
  h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
  h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
  burn = 0,  # number of burn in samples before saving posterior samples
  K = 26 # number of factors
)

MegaLMM_state = setup_model_MegaLMM(morph,
                                    ~Sexo*Poblacion + (1|Linea) + (Sexo|Linea),
                                    data = classifier,
                                    run_parameters=run_parameters,
                                    run_ID = 'MegaLMM_example')

priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.5,   nu = 5),      # Prior variance of trait residuals after accounting for fixed effects and factors
  tot_F_var = list(V = 18/20, nu = 20),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
  Lambda_prior = list(
    sampler = sample_Lambda_prec_horseshoe, # function that implements the horseshoe-based Lambda prior described in Runcie et al 2020. See code to see requirements for this function.
    prop_0 = 0.1,    # prior guess at the number of non-zero loadings in the first and most important factor
    delta = list(shape = 3, scale = 1),    # parameters of the gamma distribution giving the expected change in proportion of non-zero loadings in each consecutive factor
    delta_iterations_factor = 100   # parameter that affects mixing of the MCMC sampler. This value is generally fine.
  ),
  h2_priors_resids_fun = function(h2s,n)  1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
)

maps = make_Missing_data_map(MegaLMM_state)

MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map)
MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)  # apply the priors
MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state) # initialize the model
MegaLMM_state = initialize_MegaLMM(MegaLMM_state) # run the initial calculations

MegaLMM_state = clear_Posterior(MegaLMM_state)
MegaLMM_state = reorder_factors(MegaLMM_state)

(n_threads = optimize_n_threads(MegaLMM_state,seq(1,RcppParallel::defaultNumThreads(),by=1),times=2))
set_MegaLMM_nthreads(n_threads$optim)
# now do sampling is smallish chunks
n_iter = 100;  # how many samples to collect at once?
for(i  in 1:5) {
  print(sprintf('Run %d',i))
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_samples iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
  
  MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
  print(MegaLMM_state) # print status of current chain
  plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
  
  # set of commands to run during burn-in period to help chain converge
  if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i < 3) {
    MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    print(MegaLMM_state$run_parameters$burn)
  }
}

MegaLMM_state$Posterior = reload_Posterior(MegaLMM_state)
U_hat = get_posterior_mean(MegaLMM_state,U_R + U_F %*% Lambda)

G_samples =      get_posterior_FUN(MegaLMM_state,t(Lambda) %*% diag(F_h2['Linea',]) %*% Lambda + diag(resid_h2['Linea',]/tot_Eta_prec[1,]))

G_samples[1,,]

G_factor_samples = get_posterior_FUN(MegaLMM_state,Lambda %*% t(Lambda))

G_factor_samples[1,1:5,1:5] %>% cov2cor()


MegaLMM_state$Posterior$F_h2[,1,] %>% melt %>% 
  ggplot(., aes(as.factor(Var2), value))+geom_boxplot()





Males_Fixed<-array(Full.fit$Sol[,grep("SexoM",colnames(Full.fit$Sol))], 
                   dim = c(dim(Full.fit$Sol)[1],
                           dim(morph)[2],
                           length(unique(classifier$Poblacion))),
                   dimnames = list(NULL,
                                   colnames(morph),
                                   unique(classifier$Poblacion)))
Males_Fixed<-aperm(Males_Fixed,c(3,2,1))
for(i in 1:dim(Males_Fixed)[2]) Males_Fixed[,i,]<- 
  Males_Fixed[,i,]*attr(morph.s,"scaled:scale")[i]

mpc12<-adply(Males_Fixed[,c(1,2),], 3, 
             function(x) data.frame(pop=unique(classifier$Poblacion), x))

fixed_m.plot<-
  data.frame(morph, classifier) %>%
  # subset(., Sexo="H") %>%
  ggplot(., aes(PC1,PC2))+
  geom_point(col="gray")+
  geom_point(aes(PC1,PC2,color=pop),mpc12, alpha=0.01)+
  coord_fixed()+
  theme_bw()+
  stat_ellipse(aes(color=pop),mpc12)
ggsave("fixed_m.pdf",fixed_m.plot,width = 5,height = 4)

Females_Fixed<-array(Full.fit$Sol[,grep("SexoH",colnames(Full.fit$Sol))], 
                     dim = c(dim(Full.fit$Sol)[1],
                             dim(morph)[2],
                             length(unique(classifier$Poblacion))),
                     dimnames = list(NULL,
                                     colnames(morph),
                                     unique(classifier$Poblacion)))
Females_Fixed<-aperm(Females_Fixed,c(3,2,1))
for(i in 1:dim(Females_Fixed)[2]) Females_Fixed[,i,]<- 
  Females_Fixed[,i,]*attr(morph.s,"scaled:scale")[i]

fpc12<-adply(Females_Fixed[,c(1,2),], 3, 
             function(x) data.frame(pop=unique(classifier$Poblacion), x))

fixed_f.plot<-
  data.frame(morph, classifier) %>%
  # subset(., Sexo="H") %>%
  ggplot(., aes(PC1,PC2))+
  geom_point(col="gray")+
  geom_point(aes(PC1,PC2,color=pop),fpc12, alpha=0.01)+
  coord_fixed()+
  theme_bw()+
  stat_ellipse(aes(color=pop),fpc12)
ggsave("fixed_f.pdf",fixed_f.plot,width = 5,height = 4)

ggsave("fixeds.pdf",
       plot_grid(fixed_f.plot, 
                 fixed_m.plot,ncol = 1,labels = "AUTO"),
       width = 5,height = 8)



full_G<-array(Full.fit$VCV,
              dim = c(dim(Full.fit$VCV)[1],
                      dim(morph)[2]*2,
                      dim(morph)[2]*2),
              dimnames = list(NULL,
                              c(paste0(colnames(morph),":F"),
                                paste0(colnames(morph),":M")),
                              c(paste0(colnames(morph),":F"),
                                paste0(colnames(morph),":M"))))
full_G<-aperm(full_G,c(2,3,1)) 

for(i in 1:dim(full_G)[3]) {
  full_G[,,i]<-full_G[,,i]*
    (rep(attr(morph.s,"scaled:scale"),2)%*%
       t(rep(attr(morph.s,"scaled:scale"),2)))
}

dim.ret
mellipc1pc2<-full_G[1:2+dim.ret,1:2+dim.ret,] %>% adply(., 3, function(x){
  ellipse::ellipse(x,npoints = 50)
})
names(mellipc1pc2)[-1]<-c("PC1", "PC2")
G_m.plot<-
  data.frame(morph, classifier) %>%
  ggplot(., aes(PC1,PC2))+
  geom_point(col="gray")+
  coord_fixed()+
  theme_bw()+
  geom_polygon(aes(PC1,PC2,group=X1),mellipc1pc2, 
               fill=NA, color="black",size=0.01)
ggsave("G_m.pdf",G_m.plot,width = 4,height = 4)

fellipc1pc2<-full_G[1:2,1:2,] %>% adply(., 3, function(x){
  ellipse::ellipse(x,npoints = 50)
})
names(fellipc1pc2)[-1]<-c("PC1", "PC2")

G_f.plot<-
  data.frame(morph, classifier) %>%
  ggplot(., aes(PC1,PC2))+
  geom_point(col="gray")+
  coord_fixed()+
  theme_bw()+
  geom_polygon(aes(PC1,PC2,group=X1),fellipc1pc2, 
               fill=NA, color="black",size=0.01)
ggsave("G_f.pdf",G_f.plot,width = 4,height = 4)



save(file="Gs_11_weakPrior.Rdata", Females_Fixed, Males_Fixed, full_G, Full.fit)