names(ftregion)
ftregion<-read.csv("newftregion.txt",sep=";",header=T)
ftregion<-read.csv("ftcomuna.txt",sep=";",header=T)
ftis<-read.csv("ftiniciosintomas.txt",sep=";",header=T)
summary(ftis)
table(ftis$fch_publicacion)

ftis2=ftis[which(ftis$fch_publicacion=="2020-10-26"),]
ftis3=ftregion[which(ftregion$id_region==4),]
ftis3
names(ftis3)
View(ftis3)
table(ftis3)
ftis3=ftis2[which(ftis2$id_region==2),]
ftis3
length(ftis3)
t1=aggregate(cant_casosconfirmados~id_semanaepidemiologica,ftis3,sum)
t1
t2=aggregate(cant_casosconfirmadoschowellw5~id_semanaepidemiologica,ftis3,sum)
t2
plot(t1[,1],t1[,2],type="l")
lines(t2[,1],t2[,2],col="red")


############################
########################
##Análisis Descriptivo


###Cuerpo
 install.packages("remotes")
remotes::install_github("jespermaag/gganatogram")
library(gganatogram)
 install.packages("dplyr")
library(dplyr)

hgMale_key %>%
   filter(organ %in% c( "lung")) %>%
   gganatogram(organism = "human", sex = "male",
               fill = "colour") +
   theme_void() + 
   coord_fixed()

hgMale_key %>%
   filter(type %in% "nervous_system") %>%
   gganatogram(organism = "human", sex = "male",
               fill = "colour", outline = FALSE) +
   theme_void() + 
   coord_fixed()

##calendario

install.packages("calendR")
library(calendR)

# Datos
datos <- rnorm(30, 15, 10)

# Crea un vector donde todos los valores son ligeramente
# inferiores que el menor valor de tus datos
dias <- rep(min(datos) - 0.05, 365)

# Rellena los días que quieras con tus datos
dias[30:59] <- datos

calendR(year = 2021,
        special.days = dias,
        low.col = "white",
        special.col = "#FF0000",
        gradient = TRUE,
        legend.pos = "bottom")
###########
night_owlish <- "https://raw.githubusercontent.com/batpigandme/night-owlish/master/rstheme/night-owlish.rstheme"
rstudioapi::addTheme(night_owlish, apply = TRUE)
##########
set.seed(2)
ftregion<-read.csv("newftregion.txt",sep=";",header=T)
names(ftregion)
data_set <- data.frame(Region = ftregion$txt_nombreregion,
   poblacion = ftregion$cant_poblacion,
                       type = sample(1:4, size = 25, replace = TRUE),
                       store = sample(paste("Store", 1:4),
                                      size = 25, replace = TRUE))

head(data_set)
####################


ftregion<-read.csv("newftregion.txt",sep=";",header=T)
plot(df$e, df$b, type = "l", col ="tomato", lwd=2, ylim = c(-5,5))
lines(df$e, df$c, col = "lightblue", lwd = 2, lty= "dashed")

ftis3=ftregion[which(ftregion$id_region==0),]
a<-ftis3$id_region
e<-ftis3$fch_confirmado
p<-ftis3$cant_poblacion
b<-ftis3$cant_casosconfirmadosdiario
b=round(b/p*100000)
c<-ftis3$cant_uci
c=round(c/p*100000)
d<-ftis3$cant_falle
d=round(d/p*100000)
x<-ftis3$txt_semanaepidemiologica
df=data.frame(a,e,b,c,d,x)

ftis3
length(ftis3)
t1=aggregate(b~x,ftis3,sum)
t2=aggregate(c~x,ftis3,sum)
t3=aggregate(d~x,ftis3,sum)

qqplot(as.Date(df[,2]),df[,3],type="l")
lines(as.Date(df[,2]),df[,4],col="red")
lines(as.Date(df[,2]),df[,5],col="purple")

lines(t2[,1],t2[,2],col="red")
lines(t1[,1],t1[,2],col="red")
######################



########Trabajo de Stéphane Ghozzi asmodee-trendbreaker-evaluation.R


library(here)
library(MASS)
library(epitrix)
library(distcrete)
library(tidyr)
library(projections)
library(caret)
 remotes::install_github("reconhub/trending@bootstrap")
 remotes::install_github("reconhub/trendbreaker@bootstrap")
library(trendbreaker)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(tools)
library(ggtext)
library(surveillance)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts=F)

# A new version of the `farringtonFlexible` function written by Michael Höhle:
# - offers the option of removing the safeguard against extrapolating trend through
#   `safeguardAgainstExtrapolation=F`
# - defines bounds on one-side p-values of the prediction interval (bug fixed in
#   version 1.19.0 of package `surveillance`)
source( 'farringtonFlexible.R')


### Global parameters ----

compute_simulations <- T
compute_detections <- T
compute_scores <- T
plot_results <- T
compute_ccgs <- T
plot_ccgs <- T
download_nhs_pathways <- F


### Scenarios and simulations ----

## Parameters

# Paths
data_relative_path <- 'data'
dir.create(here(data_relative_path), showWarnings=F, recursive=T)
img_relative_path <- 'img'
dir.create(here(img_relative_path), showWarnings=F, recursive=T)
ccgs_relative_path <- paste0(img_relative_path, '/ccgs')
dir.create(here(ccgs_relative_path), showWarnings=F, recursive=T)

# Scenarios = overall simulation parameters
overall_params <- list(
   sim_methods = c('project'), # c('project', 'pointsource')
   n_sim_steps = 42L, # number of simulation steps
   n_replicates = 30L, # number of simulation per scenario and simulation method
   n_max_change = 1L, # maximum number of change points
   report_prob = 0.1, # probability of an infection being reported and thus observed
   d_min_period = 7L, # min duration of a period
   d_observation_period = 12L, # Period of observation, during which the algorithms are
   # evaluated.
   period_trends = c('strong_upward'=4, 'mild_upward'=3, 'constant'=2, 'downward'=1), # types of
   # trends within a period with relative upward strength
   change_within_observation = T, # Restrict scenarios to those with either no change in trend or at
   # least one change within the period of observation (the first step of the last period lies
   # between the first step of observation plus 3 days, `n_sim_steps-d_observation_period+4` and
   # `n_sim_steps-d_min_period+1`).
   buffer_steps = 3, # minimum number of days after start of observation after which last period
   # can start if `change_within_observation` is `TRUE` (this is to allow at least a few points
   # of different classes when `buffer_steps >= 1`).
   initial_levels = c('high', 'medium', 'low'), # initial incidence level
   classes = c('decrease','increase','normal'),
   select_interesting_scenarios = T, # Apply detection only to scenarios selected with
   # `interesting_scenarios`(or if FALSE to all scenarios)
   interesting_scenarios = list(
      steady_state=list(
         sim_method=c('project'), initial_level='medium', n_periods=1L, trends=c('constant')
      ),
      lockdown=list(
         sim_method=c('project'), initial_level='medium', n_periods=2L,
         trends=c('mild_upward','downward')
      ),
      relapse=list(
         sim_method=c('project'), initial_level='high', n_periods=2L,
         trends=c('downward','mild_upward')
      ),
      flareup=list(
         sim_method=c('project'), initial_level='low', n_periods=2L,
         trends=c('constant','strong_upward')
      )
   ), # Define scenarios which are interesting. `trends` should have length `n_periods`.
   # `interesting_scenarios` can be set to NULL if no scenario is especially interesting.
   # Even if `select_interesting_scenarios` is FALSE, `interesting_scenarios` is taken into
   # account to annotate the scenarios generated.
   detect_algos = c('ASMODEE_manual', 'ASMODEE_optimal', 'modified_Farrington'), # , 'NegBin'),
   detect_alpha = sort(c(seq(0, 0.04, by=0.01), seq(0.05, 1, by=0.05))), # alpha values on which to
   # compute scores
   alpha_plot_pod = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), # values of alpha for the
   # probability-of-detection plot, have to be in `detect_alpha`
   score_types = c('sensitivity', 'specificity', 'precision', 'f1', 'ba', 'timeliness', 'pod'),
   # "ba" = balanced accuracy, "pod" = probability of detection
   alpha_opt_type = 'pod' # Which score shall be used to set the optimal alpha? Can be either one of
   # "ba", "f1", "pod" and "sensitivity": alpha is set to maximize the score or for "sensitivity"
   # for it to as close as possible to 0.9.
)
saveRDS(overall_params, here(data_relative_path, 'overall_params.rds'))

# Parameters for `projections::project()`
project_params <- list(
   R_means = c(mild_upward=1.3, strong_upward=2.5, constant=1, downward=0.8),
   R_sd = 0.1,
   init_duration = 14L,
   init_incidence = list(high=10000L, medium=1000L, low=100L),
   si_mean = 4.7,
   si_sd = 2.9
)
saveRDS(project_params, here(data_relative_path, 'project_params.rds'))

si_param = epitrix::gamma_mucv2shapescale(4.7, 2.9/4.7)
si_distribution <- distcrete::distcrete('gamma', interval=1, shape=si_param$shape,
                                        scale=si_param$scale, w=0.5)

ProjectSimulations <- function(params) {
   inc <- incidence::incidence(rep(seq_len(params$duration_ini), each=params$n_ini))
   project(inc,
           R = params$R,
           n_sim = params$n_replicates,
           n_days = params$duration,
           time_change = params$time_change,
           si = params$si)
}

# ASMODEE parameters
asmodee_params <- list(
   k_manual = overall_params$d_observation_period,
   k_optimal_max = overall_params$d_observation_period,
   method = evaluate_aic, # evaluate_resampling,
   models = list(
      poisson_constant = glm_model(count ~ 1, family='poisson'),
      # regression = lm_model(count ~ date),
      poisson_time = glm_model(count ~ date, family = "poisson"),
      negbin_time = glm_nb_model(count ~ date)
   )
)

# modified Farrington parameters
ff_control <- list(
   b=floor((overall_params$n_sim_steps-overall_params$d_observation_period-3)/7),
   w=3,
   limit54=c(0,1),
   noPeriods=7,
   pastWeeksNotIncluded=7,
   glmWarnings=F,
   thresholdMethod='nbPlugin',
   pThresholdTrend=1,
   safeguardAgainstExtrapolation=F # new option of `farringtonFlexible` from the modified
   # farringtonFlexible.R script
)

# Consistency checks
if (overall_params$n_sim_steps < (overall_params$n_max_change+1)*overall_params$d_min_period) {
   stop(paste0(
      'Not enough simulation steps (', overall_params$n_sim_steps,
      ') given the maximum number of periods (', overall_params$n_max_change+1,
      ') and the minimum period duration (',
      overall_params$d_min_period, ').')
   )
}

if (overall_params$change_within_observation &
    overall_params$d_observation_period < overall_params$buffer_steps +
    overall_params$d_min_period) {
   stop(paste0(
      'Inconsistent choice of parameters: at least one change in the period of observation ',
      '(`overall_params$change_within_observation` is `TRUE`) but the duration of observation (',
      overall_params$d_observation_period, ') is shorter than the minimum period duration (',
      overall_params$d_min_period, ') plus the buffer (', overall_params$buffer_steps, ').',
      'Please revise.')
   )
}

if (overall_params$select_interesting_scenarios & is.null(overall_params$interesting_scenarios)) {
   stop('You selected to apply detection only to interesting scenarios but didn\'t provide any.')
}

# Colors
vega_standard_palette <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
                           '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')
major_minor_colors <- c(major=vega_standard_palette[1], minor=vega_standard_palette[2])
score_colors <- colorRampPalette(brewer.pal(9,'Pastel1'))(length(overall_params$score_types))
names(score_colors) <- overall_params$score_types
algo_colors <- colorRampPalette(brewer.pal(12,'Set3'))(length(overall_params$detect_algos))
names(algo_colors) <- overall_params$detect_algos

## Generate scenarios

# N.B. the different timings of trend change (different start dates for a given period with a
# given trend) are considered part of the same scenario, i.e. they have the same `id_scenario`.

# Trend combinations
trend_combinations <- lapply(
   1:(overall_params$n_max_change+1),
   function (np) {
      x <- expand.grid(
         as_tibble(
            replicate(np, names(overall_params$period_trends)),
            .name_repair = 'minimal'
         ),
         stringsAsFactors=F
      ) %>%
         as_tibble() %>%
         bind_cols(
            as_tibble(
               t(replicate(overall_params$n_max_change+1-np, as.character(NA))),
               .name_repair = 'minimal'
            )
         )
      names(x) <- paste0('trend_period_', 1:(overall_params$n_max_change+1))
      x
   }
)
trend_combinations <- bind_rows(trend_combinations)

# Remove combinations with same trends before and after a change in trend, as these cases are
# already covered in scenarios with one less trend change.
 trend_cols <- names(trend_combinations)[names(trend_combinations)!='n_periods']
if (ncol(trend_combinations)>=2) {
   id_remove_trendcombination <- c()
   for (j in 1:(ncol(trend_combinations)-1)) {
      id_remove_trendcombination <- c(
         id_remove_trendcombination,
         which(trend_combinations[,j+1] == trend_combinations[,j])
      )
   }
   trend_combinations <- trend_combinations[-id_remove_trendcombination,]
}
trend_combinations$n_periods <- sapply(
   1:nrow(trend_combinations),
   function (i) as.integer(overall_params$n_max_change+1-sum(is.na(trend_combinations[i,])))
)

# Possible positions of trend changes, i.e. possible first steps of each trend period.
last_start <- overall_params$n_sim_steps - overall_params$d_min_period + 1
period_starts <- lapply(
   1:(overall_params$n_max_change),
   function (nc) {
      x <- expand.grid(
         as_tibble(
            replicate(nc, 1:last_start),
            .name_repair = 'minimal'
         ),
         stringsAsFactors=F
      ) %>%
         as_tibble() %>%
         bind_cols(
            as_tibble(
               t(replicate(overall_params$n_max_change-nc, as.integer(NA))),
               .name_repair = 'minimal'
            )
         )
      names(x) <- paste0('period_start_', 2:(overall_params$n_max_change+1))
      x <- bind_cols(tibble(period_start_1 = 1L), x)
      x
   }
)
pstart_1p <- c(1L, rep(as.integer(NA), overall_params$n_max_change))
names(pstart_1p) <- paste0('period_start_', 1:(overall_params$n_max_change+1))
period_starts <- bind_rows(pstart_1p, period_starts)

# Remove trend periods shorter than the minimum duration
if (ncol(period_starts)>=2) {
   id_remove_changeposition <- c()
   for (j in 1:(ncol(period_starts)-1)) {
      id_remove_changeposition <- c(
         id_remove_changeposition,
         which(period_starts[,j+1] - period_starts[,j] < overall_params$d_min_period)
      )
   }
   period_starts <- period_starts[-id_remove_changeposition,]
}

# If so desired, remove scenarios with at least one change trend but none in the observation period.
if (overall_params$change_within_observation) {
   id_remove_not_obs <- c()
   for (i in 1:nrow(period_starts)) {
      maxchange <- max(period_starts[i,], na.rm=T)
      np <- ncol(period_starts)-sum(is.na(period_starts[i,]))
      if (np > 1 &
          maxchange < overall_params$n_sim_steps - overall_params$d_observation_period +
          overall_params$buffer_steps + 1) {
         id_remove_not_obs <- c(id_remove_not_obs, i)
      }
   }
   period_starts <- period_starts[-id_remove_not_obs,]
}

period_starts$n_periods <- sapply(
   1:nrow(period_starts),
   function (i) as.integer(overall_params$n_max_change+1-sum(is.na(period_starts[i,])))
)
period_starts <- period_starts %>% mutate(id_start_comb=row_number())
period_starts <- period_starts %>%
   pivot_longer(cols=-c(n_periods,id_start_comb), names_to='id_period', values_to='period_start') %>%
   mutate(id_period=as.integer(gsub('period_start_', '', id_period))) %>%
   filter(!is.na(period_start))

# All combinations of possible parameter values
scenarios <- expand.grid(
   sim_method = overall_params$sim_methods,
   initial_level = overall_params$initial_levels,
   n_periods = 1:(overall_params$n_max_change+1),
   stringsAsFactors=F
) %>%
   as_tibble() %>%
   full_join(trend_combinations, by='n_periods') %>%
   mutate(id_scenario=row_number()) %>%
   pivot_longer(cols=contains('trend_period_'), names_to='id_period', values_to='trend') %>%
   mutate(id_period=as.integer(gsub('trend_period_', '', id_period))) %>%
   filter(!is.na(trend)) %>%
   full_join(period_starts, by=c('n_periods', 'id_period')) %>%
   unique() %>%
   dplyr::select(id_scenario, sim_method, initial_level, n_periods, id_period, trend, id_start_comb,
                 period_start)

# Annotate scenarios with interesting ones
scenarios$interesting <- as.character(NA)
if (!is.null(overall_params$interesting_scenarios)) {
   for (isc in names(overall_params$interesting_scenarios)) {
      sm <- overall_params$interesting_scenarios[[isc]]$sim_method
      il <- overall_params$interesting_scenarios[[isc]]$initial_level
      np <- overall_params$interesting_scenarios[[isc]]$n_periods
      trends <- overall_params$interesting_scenarios[[isc]]$trend
      
      scenarios <- scenarios %>%
         mutate(interesting=replace(
            interesting,
            sim_method %in% sm & initial_level==il & n_periods==np,
            isc
         ))
      
      for (ids in scenarios %>% filter(!is.na(interesting)) %>% pull(id_scenario) %>% unique()) {
         s_trends <- scenarios %>% filter(id_scenario==ids) %>% dplyr::select(id_period, trend) %>%
            unique() %>% arrange(id_period) %>% pull(trend)
         if (!identical(s_trends, trends)) {
            scenarios <- scenarios %>%
               mutate(interesting=replace(interesting, id_scenario==ids & interesting==isc, NA))
         }
      }
   }
}

# Look only at interesting scenarios?
if (overall_params$select_interesting_scenarios) {
   scenarios <- scenarios %>% filter(!is.na(interesting))
}

saveRDS(scenarios, here(data_relative_path, 'scenarios.rds'))

### Generate simulations ----

if (compute_simulations) {
   
   simulations <- NULL
   # t1 <- Sys.time()
   for (ids in sort(unique(scenarios$id_scenario))) {
      
      ## DEBUG
      # ids <- scenarios %>% filter(interesting=='steady_state') %>% pull(id_scenario) %>% unique()
      # ids <- scenarios %>% filter(interesting=='lockdown') %>% pull(id_scenario) %>% unique()
      # ids <- scenarios %>% filter(interesting=='relapse') %>% pull(id_scenario) %>% unique()
      # ids <- scenarios %>% filter(interesting=='flareup') %>% pull(id_scenario) %>% unique()
      
      sm <- scenarios %>% filter(id_scenario==ids) %>% pull(sim_method) %>% unique()
      istart <- scenarios %>% filter(id_scenario==ids) %>% pull(id_start_comb) %>% unique()
      il <- scenarios %>% filter(id_scenario==ids) %>% pull(initial_level) %>% unique()
      trends <- scenarios %>% filter(id_scenario==ids) %>% pull(trend) %>% unique()
      interest <- scenarios %>% filter(id_scenario==ids) %>% pull(interesting) %>% unique()
      
      if (sm=='project') {
         # Simulate infections with branching process from `projections::project()`.
         # Then to simulate under-reporting: sample from the infection count: each infection is
         # reported with probability `overall_params$report_prob`.
         
         run_offset <- 0L
         for (ist in istart) {
            ## DEBUG
            # ist <- 2
            
            change_at <- scenarios %>%
               filter(id_scenario==ids & id_start_comb==ist) %>%
               pull(period_start) %>% unique()
            
            if (length(change_at)==1) {
               change_at <- NULL
            } else {
               change_at <- change_at[-1]
            }
            
            sim_param <- list(
               n_ini=project_params$init_incidence[[il]],
               R=lapply(
                  trends,
                  function(tre)
                     rlnorm(
                        overall_params$n_replicates,
                        log(project_params$R_means[[tre]]),
                        project_params$R_sd
                     )
               ),
               time_change = change_at,
               duration_ini = project_params$init_duration,
               duration = overall_params$n_sim_steps,
               n_replicates = overall_params$n_replicates,
               si = si_distribution
            )
            project_res_list <- lapply(
               1:overall_params$n_replicates,
               function (i_rep) {
                  infect_count <- as.integer(ProjectSimulations(sim_param)[,i_rep])
                  as.integer(sapply(infect_count,
                                    function (ic) rbinom(n=1, size=ic, prob=overall_params$report_prob)))
               }
            )
            names(project_res_list) <- 1:overall_params$n_replicates
            
            project_res <- bind_cols(project_res_list) %>% as_tibble()
            project_res$interesting <- interest
            project_res$id_scenario <- ids
            project_res$sim_method <- sm
            project_res$id_start_comb <- ist
            project_res$sim_step <- 1:nrow(project_res)
            project_res <- project_res %>%
               as_tibble() %>%
               pivot_longer(cols=-c(interesting, id_scenario, sim_method, id_start_comb, sim_step),
                            names_to='replicate', values_to='count') %>%
               mutate(
                  sim_run=as.integer(replicate)+run_offset,
                  replicate=as.integer(replicate),
                  class=NA) %>%
               dplyr::select(interesting, id_scenario, sim_method, id_start_comb, sim_run, replicate,
                             sim_step, count, class) %>%
               arrange(sim_run, sim_step)
            run_offset <- run_offset + overall_params$n_replicates
            
            # Classes of observations which will serve as ground truth:
            # - if there is only one period, all observations are "normal"
            # - all observations in the period before last are "normal"
            # - if the trend stays the same, then all observations of the last period are "normal"
            # - if the trend between period before last and last changes and increases
            #   (to "upward" or from "downward" to "constant"), then the observations of
            #   the last period are "increase"
            # - respectively, if the trend changes to "downward" or from "upward" to "constant",
            #   then the observations of the last period are "decrease"
            if (length(trends)==1) {
               project_res$class <- 'normal'
            } else {
               strength_t1 <- overall_params$period_trends[[trends[length(trends)-1]]]
               strength_t2 <- overall_params$period_trends[[trends[length(trends)]]]
               if (strength_t1 == strength_t2) {
                  project_res$class <- 'normal'
               } else if (strength_t1 < strength_t2) {
                  project_res <- project_res %>%
                     mutate(class=ifelse(sim_step %in% 1:(tail(change_at, 1)-1), 'normal', 'increase'))
               } else if (strength_t1 > strength_t2) {
                  project_res <- project_res %>%
                     mutate(class=ifelse(sim_step %in% 1:(tail(change_at, 1)-1), 'normal', 'decrease'))
               } else {
                  stop('Something\'s wrong with the trends: "', paste0(trends, collapse='", "'), '".')
               }
            }
            if (!all(project_res$class %in% overall_params$classes)) {
               stop('Some the attributed classes "', paste(sort(unique(project_res$class)), collapse=", "),
                    '" are not in the expected classes "', paste(overall_params$classes, collapse=", "), '".')
            }
            
            simulations <- simulations %>% bind_rows(project_res)
         }
         
      } else {
         
         stop('Don\'t know simulation method "', sm, '"')
         
      }
      
   }
   
   # t2 <- Sys.time()
   # print('Time elapsed for generating simulations:')
   # print(t2-t1)
   
   saveRDS(simulations, here(data_relative_path, 'simulations.rds'))
   
   # In the flare-up scenario, count how many time series actually go up in the
   # second period, according to 3 criteria: mean increases, median increases, or
   # the 25th percentile of the second period is larger or equal to the 75th of
   # the first.
   sim_count_diff_flareup <- simulations %>%
      filter(interesting == "flareup") %>%
      dplyr::select(sim_run, count, class) %>%
      group_by(sim_run, class) %>%
      summarize(
         mean_count = mean(count),
         med_count = median(count),
         q25_count = quantile(count, probs = 0.25),
         q75_count = quantile(count, probs = 0.75)
      ) %>%
      mutate(ref_count = case_when(
         class == "increase" ~ q25_count,
         class == "normal" ~ q75_count,
         TRUE ~ as.numeric(NA)))
   
   prop_flareup_increase <- c()
   for (countq in c("mean_count", "med_count", "ref_count")) {
      prop_flareup_increase_df <- sim_count_diff_flareup %>%
         dplyr::select(sim_run, class, all_of(countq)) %>%
         pivot_wider(names_from = class, values_from = all_of(countq)) %>%
         mutate(actual_increase = increase >= normal)
      
      prop_flareup_increase <- c(
         prop_flareup_increase,
         sum(prop_flareup_increase_df$actual_increase) /
            nrow(prop_flareup_increase_df)
      )
   }
   
   overview_relative_path <- paste0(img_relative_path, '/sim-project/overview')
   dir.create(overview_relative_path, showWarnings = FALSE, recursive = TRUE)
   prop_flareup_increase_file <- file(here(overview_relative_path,
                                           "prop_flareup_increase.txt"))
   open(prop_flareup_increase_file, open = "w")
   cat("In the flare-up scenario, proprtion of time-series that are actually",
       "going up in the second period, according to three criteria: mean increases,",
       "median increases, or the 25th percentile of the second period is larger or",
       "equal to the 75th of the first.", "\n",
       "Rt(normal) = ", project_params$R_means[[
          overall_params$interesting_scenarios$flareup$trends[1]
       ]], "\n",
       "Rt(increase) = ", project_params$R_means[[
          overall_params$interesting_scenarios$flareup$trends[2]
       ]], "\n",
       "initial incidence = ", project_params$init_incidence[[
          overall_params$interesting_scenarios$flareup$initial_level
       ]], "\n\n",
       "mean: ", signif(prop_flareup_increase[1], digits = 4), "\n",
       "median: ", signif(prop_flareup_increase[2], digits = 4), "\n",
       "percentiles: ", signif(prop_flareup_increase[3], digits = 4), "\n",
       sep = "",
       file = prop_flareup_increase_file)
   close(prop_flareup_increase_file)
   
} else {
   
   simulations <- readRDS(here(data_relative_path, 'simulations.rds'))
   
}

### Apply detection algorithms ----
#
# TODO: For NegBin and farringtonFlexible, don't loop over and retrain for each `alpha` but
#       train once and apply cut-offs on quantiles computed from the Negative Binoomial
#       distribution.

ClassifyCountCI <- function(cnt, ci) {
   # Classify observed counts `cnt` as "decrease", "normal", "increase" if they are below, within or
   # above the confidence interval `ci`.
   
   classification <- c()
   for (i in 1:length(cnt)) {
      if (is.na(ci[['lowerbound']][i]) | is.na(ci[['upperbound']][i])) {
         clas <- NA
      } else if (cnt[i] < ci[['lowerbound']][i]) {
         clas <- 'decrease'
      } else if (cnt[i] > ci[['upperbound']][i]) {
         clas <- 'increase'
      } else {
         clas <- 'normal'
      }
      classification <- c(classification, clas)
   }
   return(classification)
}

# All combinations of simulations and detection parameters leading to as many detection runs.
detection_comb <- simulations %>%
   dplyr::select(interesting, id_scenario, sim_method, id_start_comb, sim_run) %>%
   unique() %>%
   mutate(id_detect_comb=row_number()) %>%
   dplyr::select(id_detect_comb, everything())

detectmeth_alpha_comb <- expand.grid(
   id_detect_comb = detection_comb$id_detect_comb,
   detect_method = overall_params$detect_algos,
   alpha = overall_params$detect_alpha,
   stringsAsFactors=F
) %>% as_tibble()

detection_comb <- detection_comb %>%
   full_join(detectmeth_alpha_comb, by='id_detect_comb') %>%
   mutate(id_detect_comb=row_number())

if (compute_detections) {
   
   detections <- NULL
   asmodee_k <- NULL
   scenario_seen <- c()
   detect_examples <- list()
    t3 <- Sys.time()
   for (idc in detection_comb$id_detect_comb) {
      
      ## DEBUG
       idc <- detection_comb %>% filter(interesting=='relapse' & sim_method=='project' & sim_run==1 & detect_method=='ASMODEE_optimal' & alpha==0.05) %>% pull(id_detect_comb)
       idc <- detection_comb %>% filter(interesting=='relapse' & sim_method=='project' & sim_run==1 & detect_method=='NegBin' & alpha==0.05) %>% pull(id_detect_comb)
       idc <- detection_comb %>% filter(interesting=='lockdown' & sim_method=='project' & sim_run==1 & detect_method=='modified_Farrington' & alpha==0.05) %>% pull(id_detect_comb)
      
       if (round(100000*idc/nrow(detection_comb))/1000 == round(100*idc/nrow(detection_comb))) {
         cat('detect: ', idc, ' / ', round(100*idc/nrow(detection_comb)), '% // ', sep='')
         print(Sys.time()-t3)
       }
      
      ids <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(id_scenario)
      sm <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(sim_method)
      sr <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(sim_run)
      dm <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(detect_method)
      al <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(alpha)
      
      if (!ids %in% scenario_seen) {
         scenario_seen <- c(scenario_seen, ids)
         detect_examples[[ids]] <- list()
         detect_examples[[ids]][[dm]] <- list()
      }
      
      count_timeseries <- simulations %>%
         filter(id_scenario==ids & sim_run==sr) %>%
         dplyr::select(sim_step, count)
      
      if(grepl('ASMODEE', dm)) {
         if (grepl('manual', dm)) {
            asmodee_res <- asmodee(
               count_timeseries %>% rename(date=sim_step),
               models = asmodee_params$models,
               alpha = al,
               fixed_k = asmodee_params$k_manual,
               method = asmodee_params$method,
               uncertain = FALSE,
               simulate_pi = TRUE
            )
         } else if (grepl('optimal', dm)) {
            asmodee_res <- asmodee(
               count_timeseries %>% rename(date=sim_step),
               models = asmodee_params$models,
               alpha = al,
               max_k = asmodee_params$k_optimal_max,
               method = asmodee_params$method,
               uncertain = FALSE,
               simulate_pi = TRUE
            )
         } else {
            stop('Wrong ASMODEE method "', dm,
                 '". It should be either with manual k ("ASMODEE_manual") ',
                 'or optimal k ("ASMODEE_optimal").')
         }
         
         asmodee_k <- asmodee_k %>%
            bind_rows(tibble(id_detect_comb=idc, detect_method=dm, alpha=al, k=asmodee_res$k))
         
         if (sr==1 & sm=='project') {
            detect_examples[[ids]][[dm]][[as.character(al)]] <- asmodee_res
         }
         
         detection_res <- asmodee_res$results %>%
            as_tibble() %>%
            dplyr::select(date, classification) %>%
            rename(sim_step = date)
         
      } else if (dm=='NegBin') {
         
         negbin_trainset <- count_timeseries %>%
            filter(sim_step < overall_params$n_sim_steps - overall_params$d_observation_period + 1) %>%
            dplyr::select(sim_step, count)
         negbin_model <- glm.nb(count~1, link='log', data=negbin_trainset)
         if (al==0) {
            negbin_ci <- c(0,Inf)
         } else {
            nb_size <- negbin_model$theta
            nb_mu <- exp(negbin_model$coefficients[['(Intercept)']])
            negbin_ci <- list(
               lowerbound = rep(qnbinom(al/2, size=nb_size, mu=nb_mu), overall_params$n_sim_steps),
               upperbound = rep(qnbinom(1-al/2, size=nb_size, mu=nb_mu), overall_params$n_sim_steps)
            )
         }
         detection_res <- count_timeseries %>%
            mutate(classification = ClassifyCountCI(count, negbin_ci)) %>%
            dplyr::select(sim_step, classification)
         
         if (sr==1 & sm=='project') {
            negbin_res <- count_timeseries %>%
               full_join(detection_res, by='sim_step') %>%
               mutate(ci_lb=negbin_ci[['lowerbound']], ci_ub=negbin_ci[['upperbound']])
            detect_examples[[ids]][[dm]][[as.character(al)]] <- negbin_res
         }
         
      } else if (dm=='modified_Farrington') {
         
         sts <- sts(count_timeseries$count, frequency=7)
         ff_out <- farringtonFlexible(sts, control=append(list(alpha=al/2), ff_control))
         ff_size <- ff_out@control$mu0Vector/(ff_out@control$phiVector-1)
         ff_mu <- ff_out@control$mu0Vector
         ff_ci <- list(
            lowerbound = c(
               rep(NA, overall_params$n_sim_steps-length(ff_out@epoch)),
               qnbinom(al/2, size=ff_size, mu=ff_mu)
            ),
            upperbound = c(
               rep(NA, overall_params$n_sim_steps-length(ff_out@epoch)),
               qnbinom(1-al/2, size=ff_size, mu=ff_mu)
            )
         )
         detection_res <- count_timeseries %>%
            mutate(classification = ClassifyCountCI(count, ff_ci)) %>%
            dplyr::select(sim_step, classification)
         
         if (sr==1 & sm=='project') {
            ff_res <- count_timeseries %>%
               full_join(detection_res, by='sim_step') %>%
               mutate(ci_lb=ff_ci[['lowerbound']], ci_ub=ff_ci[['upperbound']])
            detect_examples[[ids]][[dm]][[as.character(al)]] <- ff_res
         }
         
      } else {
         
         stop('Don\'t know detection algorithm "', dm, '".')
         
      }
      
      detection_res <- detection_res %>%
         mutate(classification = as.character(classification)) %>%
         filter(sim_step >= overall_params$n_sim_steps - overall_params$d_observation_period + 1) %>%
         mutate(id_detect_comb=idc, id_scenario=ids, sim_method=sm, sim_run=sr,
                detect_method=dm, alpha=al) %>%
         left_join(simulations, by=c('id_scenario', 'sim_method', 'sim_run', 'sim_step')) %>%
         dplyr::select(id_detect_comb, interesting, id_scenario, sim_method, sim_run, detect_method, alpha,
                       sim_step, count, classification, class)
      
      detections <- detections %>% bind_rows(detection_res)
      
   }
   
    t4 <- Sys.time()
    print('Time elapsed for applying detection algorithms:')
    print(t4-t3)
   
   saveRDS(asmodee_k, here(data_relative_path, 'asmodee_k.rds'))
   saveRDS(detections, here(data_relative_path, 'detections.rds'))
   saveRDS(detect_examples, here(data_relative_path, 'detect_examples.rds'))
   
} else {
   
   asmodee_k <- readRDS(here(data_relative_path, 'asmodee_k.rds'))
   detections <- readRDS(here(data_relative_path, 'detections.rds'))
   detect_examples <- readRDS(here(data_relative_path, 'detect_examples.rds'))
   
}

### Compute scores ----
#
# Periods before the last are considered "normal", the last belongs to the class `period trend`.
# The evaluation is done on the last `d_max_last_period` simulation steps.

scores_considered <- expand.grid(
   score = overall_params$score_types[overall_params$score_types!='pod'],
   extent = c('micro','macro','weighted', paste('Class:', overall_params$classes)),
   stringsAsFactors=F
) %>%
   as_tibble() %>%
   filter(!(score =='timeliness' & extent %in% c('micro','macro','weighted')))

dictionary_scores <- list(
   sensitivity='Sensitivity',
   specificity='Specificity',
   precision='Pos Pred Value',
   f1 = 'f1',
   ba = 'ba',
   timeliness = 'timeliness'
)

if (compute_scores) {
   
   test_results <- NULL
    t5 <- Sys.time()
   for (idc in detection_comb$id_detect_comb) {
      
      ## DEBUG
       idc <- detection_comb %>% filter(interesting=='relapse' & sim_method=='project' & sim_run==1 & detect_method=='ASMODEE_optimal' & alpha==1) %>% pull(id_detect_comb)
       idc <- detection_comb %>% filter(interesting=='steady_state' & sim_method=='project' & sim_run==1 & detect_method=='ASMODEE_optimal' & alpha==0.4) %>% pull(id_detect_comb)
       idc <- detection_comb %>% filter(interesting=='lockdown' & sim_method=='project' & sim_run==60 & detect_method=='ASMODEE_optimal' & alpha==0.05) %>% pull(id_detect_comb)
      
       cat('sim: ', idc, ' / ', round(100*idc/nrow(detection_comb)), '%\n', sep='')
       print(Sys.time()-t5)
      
      ids <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(id_scenario)
      sm <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(sim_method)
      sr <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(sim_run)
      dm <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(detect_method)
      al <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(alpha)
      istart <- detection_comb %>% filter(id_detect_comb==idc) %>% pull(id_start_comb)
      last_start <- scenarios %>%
         filter(id_scenario==ids & id_start_comb==istart) %>%
         filter(id_period==max(id_period)) %>%
         pull(period_start)
      
      test_detection <- detections %>% filter(id_detect_comb==idc) %>% dplyr::select(classification, class)
      
      # Classification scores
      cm <- confusionMatrix(
         data=factor(test_detection$classification, levels=overall_params$classes),
         reference=factor(test_detection$class, levels=overall_params$classes)
      )
      
      scores <- cm[['byClass']]
      scores[is.nan(scores)] <- 0 # default in sklearn.metrics.precision_recall_fscore_support
      
      tp_micro <- sum(sapply(overall_params$classes, function(cl1) cm$table[cl1,cl1]), na.rm=T)
      fp_micro <- sum(sapply(overall_params$classes,
                             function(cl1) {
                                sum(sapply(overall_params$classes[overall_params$classes!=cl1],
                                           function (cl2) cm$table[cl1,cl2]), na.rm=T)
                             }), na.rm=T)
      tn_micro <- sum(sapply(overall_params$classes,
                             function(cl1) {
                                sum(sapply(overall_params$classes[overall_params$classes!=cl1],
                                           function (cl2) cm$table[cl2,cl2]), na.rm=T)
                             }), na.rm=T)
      fn_micro <- sum(sapply(overall_params$classes,
                             function(cl1) {
                                sum(sapply(overall_params$classes[overall_params$classes!=cl1],
                                           function (cl2) cm$table[cl2,cl1]), na.rm=T)
                             }), na.rm=T)
      
      class_relative_weights <- sapply(1:ncol(cm$tabl), function (j) sum(cm$table[,j]))
      class_relative_weights <- class_relative_weights/sum(class_relative_weights)
      
      for (sc in 1:nrow(scores_considered)) {
         
         ## DEBUG
          sc <- which(scores_considered$score=='timeliness' & scores_considered$extent=='Class: normal')
         
         extent <- scores_considered[sc,'extent'][[1]]
         score_type <- scores_considered[sc,'score'][[1]]
         
         if (score_type == 'timeliness') {
            
            # Timeliness
            start_observation <- overall_params$n_sim_steps - overall_params$d_observation_period + 1
            ref_period_start <- ifelse(
               last_start >= start_observation,
               last_start,
               start_observation
            )
            
            period_class <- tail(test_detection$class, 1)
            if (any(test_detection$classification %in% period_class)) {
               first_detection <- start_observation +
                  min(which(test_detection$classification == period_class)) - 1
               score <- ifelse(
                  first_detection < ref_period_start,
                  NA,
                  first_detection - ref_period_start
               )
            } else {
               score <- NA
            }
            
         } else if (score_type %in% c('sensitivity', 'specificity', 'precision')) {
            
            # Classification scores
            if(extent=='micro') {
               score <- switch(score_type,
                               sensitivity = tp_micro/(tp_micro+fn_micro),
                               specificity = tn_micro/(tn_micro+fp_micro),
                               precision = tp_micro/(tp_micro+fp_micro)
               )
            } else if(extent=='macro') {
               score <- mean(scores[,dictionary_scores[[score_type]]], na.rm=T)
            } else if(extent=='weighted') {
               score <- mean(scores[,dictionary_scores[[score_type]]]*class_relative_weights, na.rm=T)
            } else { # class specific score
               score <- scores[extent,dictionary_scores[[score_type]]]
            }
            
         }
         # When no class has a defined score, i.e. score = NA, then the mean score is NaN...
         # Replace the value by NA for better handling.
         score <- ifelse(is.nan(score), NA, score)
         
         test_scores <- tibble(
            id_detect_comb = idc,
            id_scenario = ids,
            sim_method = sm,
            sim_run = sr,
            detect_method = dm,
            alpha = al,
            score_type = score_type,
            extent = extent,
            score_value = score
         )
         
         test_results <- test_results %>% bind_rows(test_scores)
         
      }
   }
   
   # F1 score = 2*precision*sensitivity/(precision+sensitivity)
   # Balanced accuracy = (sensitivity+specificity)/2 (following binary definition... definition
   #   for multi-class is the average of sensitivity, so sensitivity/macro here)
   test_results <- test_results %>%
      pivot_wider(names_from=score_type, values_from=score_value) %>%
      mutate(
         f1 = replace(2*precision*sensitivity/(precision+sensitivity), precision==0 & sensitivity==0, 0),
         ba = (sensitivity+specificity)/2) %>%
      pivot_longer(cols=c(sensitivity, specificity, precision, f1, ba, timeliness),
                   names_to='score_type', values_to='score_value')
   
   # t6 <- Sys.time()
   # print('Time elapsed for computing scores:')
   # print(t6-t5)
   
   saveRDS(test_results, here(data_relative_path, 'test_results.rds'))
   
} else {
   
   test_results <- readRDS(here(data_relative_path, 'test_results.rds'))
   
}

### Plot results ----

# For each scenario, a figure of combined plots:
# - all simulation runs in one plot
# Then one column for each algorithm:
# - one example of detection (on run 1, project and alpha = 0.05)
# - scores (average over runs) vs. alpha
# - ROC curve (average over runs) with AUC
# - precision-recall curve (average over runs)
# - score distributions for alpha = 0.05

if (plot_results) {
   
   # t7 <- Sys.time()
   plots_overview_list <- list()
   plots_combined_eval_list <- list()
   timeliness_cum_distrib_df_all <-NULL
   scores_alpha_df_all <-NULL
   optimal_alphas_all <-NULL
   for (sm in sort(unique(detection_comb$sim_method))) {
      
      overview_relative_path <- paste0(img_relative_path, '/sim-', sm, '/overview')
      dir.create(here(overview_relative_path), showWarnings=F, recursive=T)
      file_suffix_overview <- paste0('_sim-', sm)
      
      plots_overview_list[[sm]] <- list()
      plots_combined_eval_list[[sm]] <- list()
      plots_combined_eval_list[[sm]][['timeliness']] <- list()
      plots_combined_eval_list[[sm]][['score_distrib']] <- list()
      
      for (ids in sort(unique(detection_comb$id_scenario))) {
         
         interest <- scenarios %>% filter(id_scenario==ids) %>% pull(interesting) %>% unique()
         
         compare_algos_relative_path <- paste0(img_relative_path,
                                               '/sim-', sm, '/scenario-', ids, ifelse(is.na(interest), '', paste0('-', interest)))
         dir.create(here(compare_algos_relative_path), showWarnings=F, recursive=T)
         file_suffix_compare_algos <- paste0('_sim-', sm, '_scenario-', ids,
                                             ifelse(is.na(interest), '', paste0('-', interest)))
         
         ## DEBUG
          ids <- scenarios %>% filter(interesting=='steady_state') %>% pull(id_scenario) %>% unique()
          ids <- scenarios %>% filter(interesting=='lockdown') %>% pull(id_scenario) %>% unique()
          ids <- scenarios %>% filter(interesting=='relapse') %>% pull(id_scenario) %>% unique()
          sm <- 'project'
          dm <- 'ASMODEE_optimal'
          dm <- 'NegBin'
         
         classes <- simulations %>% filter(id_scenario==ids) %>% pull(class) %>% unique()
         class_last_period <- simulations %>%
            filter(id_scenario==ids & sim_step==max(sim_step)) %>%
            pull(class) %>%
            unique()
         p_starts <- scenarios %>%
            filter(id_scenario==ids) %>%
            pull(period_start) %>%
            unique()
         
         # Simulation runs
         sim_run_df <- simulations %>%
            filter(id_scenario==ids & sim_method==sm) %>%
            dplyr::select(sim_run, sim_step, count)
         sim_run_plot <- ggplot(sim_run_df, aes(sim_step, count, color=as.factor(sim_run))) +
            geom_line(size=0.3, alpha=0.2) +
            scale_color_manual(values=rep(major_minor_colors[['major']],
                                          length(unique(sim_run_df$sim_run)))) +
            geom_vline(
               xintercept=overall_params$n_sim_steps - overall_params$d_observation_period + 1 - 0.6,
               size=0.75, color='grey80') +
            geom_vline(xintercept=p_starts - 0.6, size=0.5, linetype='dotted',
                       color=major_minor_colors[['minor']]) +
            ggtitle(
               paste('scenario: ', gsub('_', ' ', ifelse(is.na(interest), ids,
                                                         ifelse(interest=='flareup', 'flare-up', interest)))),
               subtitle='simulated time series') +
            xlab('day') +
            theme_bw() +
            theme(legend.position='none')
         plots_overview_list[[sm]][[paste0(as.character(ids),'-sims')]] <- sim_run_plot
         
         # Average scores as a function of alpha
         timeliness_df <- test_results %>%
            filter(id_scenario==ids & sim_method==sm &
                      score_type=='timeliness' &
                      (score_value<=overall_params$d_min_period-1 | is.na(score_value)) &
                      extent==paste('Class:', class_last_period) &
                      detect_method %in% overall_params$detect_algos)
         n_sim_run <- simulations %>% filter(id_scenario==ids & sim_method==sm) %>%
            pull(sim_run) %>% unique() %>% length()
         timeliness_cum_distrib_df <- timeliness_df %>%
            dplyr::select(id_scenario, sim_method, detect_method, score_value, alpha) %>%
            count(id_scenario, sim_method, detect_method, score_value, alpha) %>%
            complete(id_scenario, sim_method, detect_method, score_value, alpha, fill=list(n=0)) %>%
            filter(!is.na(score_value)) %>%
            mutate(freq=n/n_sim_run) %>%
            arrange(detect_method, alpha, score_value) %>%
            group_by(id_scenario, sim_method, detect_method, alpha) %>%
            mutate(cum_freq=cumsum(freq)) %>%
            ungroup()
         scores_alpha_df <- test_results %>%
            filter(id_scenario==ids & sim_method==sm & score_type!='timeliness' &
                      extent %in% c(paste('Class:', classes), 'micro', 'macro', 'weighted')) %>%
            group_by(id_scenario, sim_method, detect_method, score_type, extent, alpha) %>%
            summarise(score_value=mean(score_value), .groups='drop') %>%
            bind_rows(
               timeliness_cum_distrib_df %>%
                  filter(score_value==max(score_value)) %>%
                  mutate(score_value=cum_freq, extent=paste('Class:', class_last_period),
                         score_type='pod') %>%
                  dplyr::select(id_scenario, sim_method, detect_method, score_type, extent, alpha, score_value)
            )
         timeliness_cum_distrib_df_all <- timeliness_cum_distrib_df_all %>%
            bind_rows(timeliness_cum_distrib_df)
         scores_alpha_df_all <- scores_alpha_df_all %>% bind_rows(scores_alpha_df)
         
         # Optimal alpha's for each scenario and simulation method:
         #ç alpha(pod) = maximize probability of detection within `overall_params$d_min_period` days
         # alpha(ba) = maximize balanced accuracy of class of interest
         # alpha(f1) = maximize F1 of class of interest
         # alpha(sensitivity) = sensitivity of class of interest closest to 0.9
         # When there are many optimal alpha's, the smallest value is taken so as to be conservative
         # and avoid subsequent false increase or decrease alarms.
          alpha=0
          #and alpha=1 are excluded as potential optima as they have no practical value.
         optimal_alphas <- NULL
         for (sty in c('pod', 'ba', 'f1', 'sensitivity')) {
            ## DEBUG
             sty <- 'pod'
             sty <- 'sensitivity'
            alpha_opt <- scores_alpha_df %>% filter(score_type==sty, !alpha %in% c(0,1))
            if (sty=='sensitivity') {
               alpha_opt <- alpha_opt %>% mutate(score_value_opt=-abs(score_value-0.9))
            } else {
               alpha_opt <- alpha_opt %>% mutate(score_value_opt=score_value)
            }
            alpha_opt <- alpha_opt %>%
               group_by(detect_method, score_type, extent) %>%
               filter(score_value==max(score_value)) %>%
               filter(alpha==min(alpha)) %>%
               ungroup() %>%
               dplyr::select(id_scenario, sim_method, detect_method, score_type, extent, alpha, score_value)
            optimal_alphas <- optimal_alphas %>% bind_rows(alpha_opt)
         }
         optimal_alphas <- optimal_alphas %>%
            arrange(id_scenario, sim_method, detect_method, score_type, extent)
         optimal_alphas_all <- optimal_alphas_all %>% bind_rows(optimal_alphas)
         
         scores_alpha_plot <- ggplot(scores_alpha_df %>%
                                        mutate(detect_method=gsub('_',' ',detect_method)) %>%
                                        rename(`score type`=score_type, score=score_value),
                                     aes(x=alpha, y=score, color=`score type`, fill=`score type`)) +
            geom_line(color='black', size=0.25) +
            geom_point(shape=21, color='black') +
            geom_point(data=optimal_alphas %>% filter(extent==paste('Class:',class_last_period) &
                                                         score_type==overall_params$alpha_opt_type) %>%
                          mutate(detect_method=gsub('_',' ',detect_method)) %>%
                          rename(`score type`=score_type, score=score_value),
                       shape=1, color='black', size=3, show.legend=F) +
            scale_color_manual(values=score_colors) +
            scale_fill_manual(values=score_colors) +
            ylim(c(0,1)) +
            ggtitle(paste('scenario: ', gsub('_', ' ', ifelse(is.na(interest), ids,
                                                              ifelse(interest=='flareup', 'flare-up', interest)))),
                    subtitle = paste0("Multiclass classification: scores restricted on ",
                                      "different periods (first column(s)) or on three classes at the ",
                                      "same time (last three columns)")) +
            facet_wrap(detect_method ~ extent, nrow=length(overall_params$detect_algos)) +
            theme_bw()
         ggsave(scores_alpha_plot, filename=here(compare_algos_relative_path, paste0('scores_alpha',
                                                                                     file_suffix_compare_algos, '.pdf')),
                height=22, width=35, units='cm')
         ggsave(scores_alpha_plot, filename=here(compare_algos_relative_path, paste0('scores_alpha',
                                                                                     file_suffix_compare_algos, '.png')),
                height=22, width=35, units='cm', dpi=150)
         
         # Detection examples
         for (dm in sort(unique(detection_comb$detect_method))) {
            
            istart <- simulations %>%
               filter(id_scenario==ids & sim_method==sm & sim_run==1) %>%
               pull(id_start_comb) %>% unique()
            pstart <-period_starts %>% filter(id_start_comb==istart) %>% pull(period_start)
            alpha_opt_dm <- optimal_alphas %>%
               filter(detect_method==dm & score_type==overall_params$alpha_opt_type &
                         extent==paste('Class:', class_last_period)) %>%
               pull(alpha)
            detect_example_plot_df <- detect_examples[[ids]][[dm]][[as.character(alpha_opt_dm)]]
            if (grepl('ASMODEE', dm)) {
               detect_example_plot <- plot(detect_example_plot_df, 'date', point_size=1, guide=F)
            } else if (dm %in% c('modified_Farrington', 'NegBin')) {
               detect_example_plot <-
                  ggplot(detect_example_plot_df, aes(x=sim_step, y=count), point_size=1, guide=F)
               if (dm=='NegBin') {
                  detect_example_plot <- detect_example_plot +
                     geom_vline(
                        xintercept=overall_params$n_sim_steps-overall_params$d_observation_period+1-0.5,
                        linetype=2
                     )
               }
               detect_example_plot <- detect_example_plot +
                  geom_ribbon(aes(ymin=ci_lb, ymax=ci_ub), alpha = 0.4, fill='#BBB67E') +
                  geom_point(aes(color=classification), size = 1 +
                                detect_example_plot_df %>% mutate(psize=ifelse(classification=='normal', 0, 1)) %>%
                                pull(psize)) +
                  geom_line(alpha=0.3, color='black') +
                  scale_color_manual(values=list(normal='#8B8B8C', increase='#CB3355',
                                                 decrease='#32AB96')) +
                  theme_bw() +
                  theme(legend.position='none')
            }
            detect_example_plot <- detect_example_plot +
               geom_vline(
                  xintercept=overall_params$n_sim_steps - overall_params$d_observation_period + 1 - 0.6,
                  size=0.75, color='grey80') +
               geom_vline(xintercept=pstart - 0.6, size=0.5, linetype='dotted',
                          color=major_minor_colors[['minor']]) +
               ggtitle('', subtitle=paste(gsub('_', ' ', dm), '/ alpha =', alpha_opt_dm)) +
               xlab('day')
            
            plots_overview_list[[sm]][[paste0(as.character(ids),'-',dm)]] <- detect_example_plot
         }
         
         # Probability of detection vs. delay to detection
         # for selected values of alpha in `overall_params$alpha_plot_pod`
         timeliness_alpha_colors <- colorRampPalette(
            c(vega_standard_palette[1], vega_standard_palette[2]))(length(overall_params$alpha_plot_pod)
            )
         timeliness_alpha_plot <- ggplot(timeliness_cum_distrib_df %>%
                                            filter(alpha %in% overall_params$alpha_plot_pod) %>%
                                            mutate(detect_method=gsub('_',' ',detect_method)),
                                         aes(x=score_value, y=cum_freq, color=as.factor(alpha))) +
            geom_line() +
            geom_point(shape=16) +
            scale_color_manual(values=timeliness_alpha_colors) +
            scale_fill_manual(values=timeliness_alpha_colors) +
            ylim(c(0,1)) +
            facet_wrap(~detect_method, ncol=length(overall_params$detect_algos)) +
            labs(x='detection interval (days)', y='probability', color='alpha') +
            ggtitle(paste('scenario: ', gsub('_', ' ', ifelse(is.na(interest), ids,
                                                              ifelse(interest=='flareup', 'flare-up', interest)))),
                    subtitle=paste0('evaluation on class "', class_last_period, '"')) +
            theme_bw()
         ggsave(timeliness_alpha_plot, filename=here(compare_algos_relative_path,
                                                     paste0('timeliness_alpha', file_suffix_compare_algos, '.pdf')),
                height=15, width=25, units='cm')
         ggsave(timeliness_alpha_plot, filename=here(compare_algos_relative_path,
                                                     paste0('timeliness_alpha', file_suffix_compare_algos, '.png')),
                height=15, width=25, units='cm', dpi=150)
         # For the overview, discard scenario "steady state" for the plot as the scores are trivial
         if (interest != 'steady_state') {
            plots_combined_eval_list[[sm]][['timeliness']][[as.character(ids)]] <- timeliness_alpha_plot
         }
         
         # "ROC curve"
         roc_curve_df <- scores_alpha_df %>%
            filter(score_type %in% c('sensitivity', 'specificity')) %>%
            pivot_wider(names_from=score_type, values_from=score_value) %>%
            arrange(alpha)
         roc_curve_plot <- ggplot(roc_curve_df %>%
                                     mutate(detect_method=gsub('_',' ',detect_method)),
                                  aes(x=1-specificity, y=sensitivity)) +
            geom_path() +
            geom_point(aes(fill=alpha), color='black', shape=21) +
            labs(x='false positive rate (1-specificity)', y='true positive rate (sensitivity)') +
            xlim(c(0,1)) +
            ylim(c(0,1)) +
            facet_wrap(detect_method ~ extent, nrow=length(overall_params$detect_algos)) +
            ggtitle(paste('scenario: ', gsub('_', ' ', ifelse(is.na(interest), ids,
                                                              ifelse(interest=='flareup', 'flare-up', interest)))),
                    subtitle = paste0("Multiclass classification: scores restricted on ",
                                      "different periods (first column(s)) or on three classes at the ",
                                      "same time (last three columns)")) +
            theme_bw()
         ggsave(roc_curve_plot, filename=here(compare_algos_relative_path,
                                              paste0('tprfpr_curve', file_suffix_compare_algos, '.pdf')),
                height=22, width=35, units='cm')
         ggsave(roc_curve_plot, filename=here(compare_algos_relative_path,
                                              paste0('tprfpr_curve', file_suffix_compare_algos, '.png')),
                height=22, width=35, units='cm', dpi=150)
         
         # Precision-recall curve (average over runs)
         precrec_curve_df <- scores_alpha_df %>%
            filter(score_type %in% c('sensitivity', 'precision')) %>%
            pivot_wider(names_from=score_type, values_from=score_value) %>%
            arrange(alpha)
         precrec_curve_plot <- ggplot(precrec_curve_df %>%
                                         mutate(detect_method=gsub('_',' ',detect_method)),
                                      aes(x=sensitivity, y=precision)) +
            geom_path() +
            geom_point(aes(fill=alpha), color='black', shape=21) +
            labs(x='recall (sensitivity)') +
            xlim(c(0,1)) +
            ylim(c(0,1)) +
            facet_wrap(detect_method ~ extent, nrow=length(overall_params$detect_algos)) +
            ggtitle(paste('scenario: ', gsub('_', ' ', ifelse(is.na(interest), ids,
                                                              ifelse(interest=='flareup', 'flare-up', interest)))),
                    subtitle = paste0("Multiclass classification: scores restricted on ",
                                      "different periods (first column(s)) or on three classes at the ",
                                      "same time (last three columns)")) +
            theme_bw()
         ggsave(precrec_curve_plot, filename=here(compare_algos_relative_path,
                                                  paste0('precrec_curve', file_suffix_compare_algos, '.pdf')),
                height=22, width=35, units='cm')
         ggsave(precrec_curve_plot, filename=here(compare_algos_relative_path,
                                                  paste0('precrec_curve', file_suffix_compare_algos, '.png')),
                height=22, width=35, units='cm', dpi=150)
         
         # Score distributions over the runs
         alpha_opt <- optimal_alphas %>%
            filter(score_type==overall_params$alpha_opt_type) %>%
            dplyr::select(detect_method, alpha)
         score_distrib_df <- test_results %>%
            filter(id_scenario==ids & sim_method==sm & score_type!='timeliness' &
                      score_type %in% overall_params$score_types &
                      detect_method %in% overall_params$detect_algos) %>%
            right_join(alpha_opt, by='detect_method', suffix=c('', '_optimal')) %>%
            filter(alpha==alpha_optimal) %>%
            dplyr::select(sim_run, score_type, extent, detect_method, score_value)
         
         for (whichextents in c('all', 'period_class')) {
            score_distrib_df_whichextents <- score_distrib_df
            if (whichextents=='period_class') {
               score_distrib_df_whichextents <- score_distrib_df_whichextents %>%
                  filter(extent==paste('Class:', class_last_period))
            }
            score_distrib_plot <- ggplot(score_distrib_df_whichextents %>%
                                            mutate(detect_method=gsub('_',' ',detect_method)),
                                         aes(x=detect_method, y=score_value, fill=score_type)) +
               geom_violin(color='grey50', trim=T, draw_quantiles=c(0.25, 0.5, 0.75),
                           alpha=0.1) +
               geom_point(data=score_distrib_df_whichextents %>%
                             group_by(score_type, extent, detect_method) %>%
                             summarize(score_value=mean(score_value)) %>%
                             ungroup() %>%
                             mutate(detect_method=gsub('_',' ',detect_method)),
                          shape=21, color='black', size=3, stroke=1) +
               scale_fill_manual(values=score_colors) +
               ylim(c(0,1)) +
               labs(x='detection method', y='score') +
               ggtitle(paste('scenario: ', gsub('_', ' ', ifelse(is.na(interest), ids,
                                                                 ifelse(interest=='flareup', 'flare-up', interest)))),
                       subtitle=paste0(
                          ifelse(whichextents=='period_class',
                                 paste0('evaluation on class "', class_last_period, '"\n'), ''),
                          paste(
                             sapply(1:nrow(alpha_opt),
                                    function (i)
                                       paste0('alpha(', gsub('_',' ',alpha_opt$detect_method[i]), ') = ',
                                              alpha_opt$alpha[i]
                                       )
                             ), collapse=', '
                          )
                       )
               ) +
               theme_bw() +
               theme(legend.position='none', axis.text.x=element_text(angle=45, hjust=1))
            
            if (whichextents=='period_class') {
               score_distrib_plot <- score_distrib_plot +
                  facet_wrap(~score_type,
                             ncol=length(
                                overall_params$score_types[!overall_params$score_types %in% c('timeliness','pod')]
                             )
                  )
               ggsave(score_distrib_plot, filename=here(compare_algos_relative_path,
                                                        paste0('score_distrib-', whichextents, file_suffix_compare_algos, '.pdf')),
                      height=12, width=25, units='cm')
               ggsave(score_distrib_plot, filename=here(compare_algos_relative_path,
                                                        paste0('score_distrib-', whichextents, file_suffix_compare_algos, '.png')),
                      height=12, width=25, units='cm', dpi=150)
               # For the overview, discard scenario "steady state" as the scores are trivial
               if (interest != 'steady_state') {
                  plots_combined_eval_list[[sm]][['score_distrib']][[as.character(ids)]] <-
                     score_distrib_plot
               }
            } else {
               score_distrib_plot <- score_distrib_plot +
                  facet_wrap(extent~score_type,
                             ncol=length(
                                overall_params$score_types[!overall_params$score_types %in% c('timeliness','pod')]
                             )
                  )
               ggsave(score_distrib_plot, filename=here(compare_algos_relative_path,
                                                        paste0('score_distrib-', whichextents, file_suffix_compare_algos, '.pdf')),
                      height=30, width=25, units='cm')
               ggsave(score_distrib_plot, filename=here(compare_algos_relative_path,
                                                        paste0('score_distrib-', whichextents, file_suffix_compare_algos, '.png')),
                      height=30, width=25, units='cm', dpi=150)
            }
         }
      }
   }
   
   # Overview plots
   for (sm in overall_params$sim_methods) {
      
      plots_overview <- arrangeGrob(grobs=plots_overview_list[[sm]], as.table=T,
                                    ncol=length(overall_params$detect_algos)+1)
      ggsave(plots_overview, filename=here(overview_relative_path,
                                           paste0('overview-timeseries', file_suffix_overview, '.pdf')),
             width=30, height=25, unit='cm')
      ggsave(plots_overview, filename=here(overview_relative_path,
                                           paste0('overview-timeseries', file_suffix_overview, '.png')),
             width=30, height=25, unit='cm', dpi=150)
      
      plots_combined_timeliness <- arrangeGrob(
         grobs=plots_combined_eval_list[[sm]][['timeliness']],
         as.table=F,
         nrow=length(plots_combined_eval_list[[sm]][['timeliness']]))
      ggsave(plots_combined_timeliness, filename=here(overview_relative_path,
                                                      paste0('overview-timeliness', file_suffix_overview, '.pdf')),
             width=20, height=20, unit='cm')
      ggsave(plots_combined_timeliness, filename=here(overview_relative_path,
                                                      paste0('overview-timeliness', file_suffix_overview, '.png')),
             width=20, height=20, unit='cm', dpi=150)
      
      plots_combined_score_distrib <- arrangeGrob(
         grobs=plots_combined_eval_list[[sm]][['score_distrib']],
         as.table=F,
         nrow=length(plots_combined_eval_list[[sm]][['score_distrib']]))
      ggsave(plots_combined_score_distrib, filename=here(overview_relative_path,
                                                         paste0('overview-score_distrib', file_suffix_overview, '.pdf')),
             width=23, height=27, unit='cm')
      ggsave(plots_combined_score_distrib, filename=here(overview_relative_path,
                                                         paste0('overview-score_distrib', file_suffix_overview, '.png')),
             width=23, height=27, unit='cm', dpi=150)
   }
   saveRDS(timeliness_cum_distrib_df_all, here(data_relative_path, 'timeliness_cum_distrib_df_all.rds'))
   saveRDS(scores_alpha_df_all, here(data_relative_path, 'scores_alpha_df_all.rds'))
   saveRDS(optimal_alphas_all, here(data_relative_path, 'optimal_alphas_all.rds'))
   
    t8 <- Sys.time()
    print('Time elapsed for plotting:')
    print(t8-t7)
   
} else {
   
   optimal_alphas_all <- readRDS(here(data_relative_path, 'optimal_alphas_all.rds'))
   
}

### Illustration on NHS pathways of English CCG's ----

CleanCCGName <- function (name_vec) {
   ccg_names <- sapply(name_vec, function (sc)
      toTitleCase(paste(strsplit(sc, '_')[[1]], collapse=' '))
   )
   ccg_names <- gsub('Nhs ', '', ccg_names)
   return(ccg_names)
}

# Reproducing and following https://github.com/thibautjombart/nhs_pathways_monitoring/blob/master/content/post/2020-05-31-analyses-ccg.Rmd
# Published here: https://covid19-nhs-pathways-asmodee.netlify.app/

# Parameters
# Duration of Leicester cluster suggested to be at least 11-24 June:
#     https://www.bbc.com/news/uk-england-leicestershire-53257835
# Lockdown was imposed on 29 June 2020 and covers Leicester City CCG" nhs_leicester_city" as well
# as small (but populous?) parts of East Leicestershire and Rutland
# "nhs_east_leicestershire_and_rutland" and West Leicestershire "nhs_west_leicestershire":
#     https://in.reuters.com/article/health-coronavirus-britain-leicester-idINL8N2E66DY
#     https://www.bbc.com/news/uk-england-leicestershire-53229371
#     https://geoportal.statistics.gov.uk/datasets/clinical-commissioning-groups-april-2020-full-clipped-boundaries-en?geometry=-1.698%2C52.559%2C-0.544%2C52.705
# However these two surrounding CCG's don't show a stark uptick in cases.
#
# We also look at Blackburn "nhs_blackburn_with_darwen". It has had increased restrictions from
# 25 July onward:
#     https://www.gov.uk/government/news/pausing-of-lockdown-easements-in-blackburn-with-darwen-and-luton
# and a lockdown from 9 August onward:
#     https://www.lancashiretelegraph.co.uk/news/18629208.rules-lockdown-blackburn-parts-e-lancs/
#
# Greater Manchester has experienced a larger outbreak at the end of July but none of the CCG's
# looked at showed marked outbreaks in the case counts ("nhs_manchester", "nhs_bury",
# "nhs_salford" and "nhs_trafford").

selected_ccgs <- c('nhs_leicester_city', 'nhs_blackburn_with_darwen')
# 'nhs_east_leicestershire_and_rutland', 'nhs_west_leicestershire', 'nhs_manchester',
# 'nhs_bury', 'nhs_salford', 'nhs_trafford'
lockdown_start <- list(
   nhs_leicester_city = as.Date('2020-06-29'),
   nhs_blackburn_with_darwen = as.Date('2020-08-09')
)

if (compute_ccgs) {
   
    t9 <- Sys.time()
   date_range_selected <- seq(as.Date('2020-06-01'), as.Date('2020-08-10'), by='day')
   alpha_nhs <- 0.05
   asmodee_conf_vec <- c('manual_7') # , 'manual_12', 'opt')
   
   # define candidate models
   models <- list(
      # regression = lm_model(count ~ day),
      poisson_constant = glm_model(count ~ 1, family = "poisson"),
      poisson_time = glm_model(count ~ day, family = "poisson"),
      negbin_time = glm_nb_model(count ~ day),
      negbin_time_weekday = glm_nb_model(count ~ day + weekday),
      negbin_time_weekday2 = glm_nb_model(count ~ day * weekday)
   )
   
   if (download_nhs_pathways) {
      # download data
      pathways <- tempfile()
      download.file(
         "https://github.com/qleclerc/nhs_pathways_report/raw/master/data/rds/pathways_latest.rds",
         pathways)
      pathways <- readRDS(pathways)
      pathways <- as_tibble(pathways) %>%
         mutate(ccg_name = sub("_ccg$", "", ccg_name))
      counts_ccg_calls <- pathways %>%
         filter(site_type %in% c("111", "999")) %>%
         group_by(ccg_name, date, day, weekday) %>%
         summarise(count = sum(count)) %>%
         complete(date, fill = list(count = 0)) %>%
         split(.$ccg_name)
      saveRDS(counts_ccg_calls, here(data_relative_path, 'counts_ccg_calls.rds'))
   } else {
      counts_ccg_calls <- readRDS(here(data_relative_path, 'counts_ccg_calls.rds'))
   }
   # Remove the CCG "null" from the list
   counts_ccg_calls[['null']] <- NULL
   
   asmodee_selected_plot <- list()
   res_ccg_calls_list <- list()
   ccg_rank_date_all <- NULL
   
   # Try three configurations for ASMODEE:
   # - manual_7: with fixed k=7, method `evaluate_aic`, a period of observation of 7 days and
     alpha=0.05
   # - manual_12: with fixed k=12, method `evaluate_aic`, a period of observation of 12 days and
   #   alpha=0.05
   # - optimal: with optimal k, the same method as in the simulations, a period of observation of
   #   `overall_params$d_observation_period`, as in the simulations and alpha taking the value
   #   that maximizes the probability of detection in the "relapse" simulation scenario above
   for (asmodee_conf in asmodee_conf_vec) {
      
      res_ccg_calls_list[[asmodee_conf]] <- list()
      asmodee_selected_plot[[asmodee_conf]] <- list()
      ccg_rank_date <- NULL
      
      for (i in 1:length(date_range_selected)) {
         
         last_date <- date_range_selected[i]
         first_date <- last_date - overall_params$n_sim_steps + 1
         
         # Keep only CCG's with enough dates to apply ASMODEE.
         ccg_enough_dates <- which(sapply(
            counts_ccg_calls,
            function (cccalls)
               min(cccalls$date) <= first_date & max(cccalls$date) >= last_date
         ))
         counts_ccg_calls_enough_dates <- counts_ccg_calls[ccg_enough_dates]
         
         ## analyses by CCG
         ## note: for all 210 CCG's this takes about 1 minute to run with AIC model selection
         ##   (method=evaluate_aic), around 16-17 min with cross validation
         ##   (method=evaluate_resampling, the default)
         if (asmodee_conf=='manual_7') {
            obs_period <- 7
            res_ccg_calls <- lapply(counts_ccg_calls_enough_dates,
                                    function (cccalls)
                                       asmodee(
                                          data = cccalls %>% filter(date >= first_date & date <= last_date),
                                          models = models,
                                          method = asmodee_params$method,
                                          fixed_k = 7,
                                          alpha = alpha_nhs,
                                          uncertain = FALSE,
                                          simulate_pi = TRUE
                                       )
            )
         } else if (asmodee_conf=='manual_12') {
            obs_period <- 12
            res_ccg_calls <- lapply(counts_ccg_calls_enough_dates,
                                    function (cccalls)
                                       asmodee(
                                          data = cccalls %>% filter(date >= first_date & date <= last_date),
                                          models = models,
                                          method = asmodee_params$method,
                                          fixed_k = 12,
                                          alpha = alpha_nhs,
                                          uncertain = FALSE,
                                          simulate_pi = TRUE
                                       )
            )
         } else if (asmodee_conf=='opt') {
            obs_period <- overall_params$d_observation_period
            res_ccg_calls <- lapply(counts_ccg_calls_enough_dates,
                                    function (cccalls)
                                       asmodee(
                                          data = cccalls %>% filter(date >= first_date & date <= last_date),
                                          models = models,
                                          max_k = asmodee_params$k_optimal_max,
                                          method = asmodee_params$method,
                                          alpha = alpha_nhs,
                                          uncertain = FALSE,
                                          simulate_pi = TRUE
                                       )
            )
         }
         res_ccg_calls_list[[asmodee_conf]][[i]] <- res_ccg_calls
         
         # Compute CCG ranks with respect to the number of "increase" outliers
         ccg_rank <- lapply(res_ccg_calls, function(e)
            tibble(
               ccg_name = unique(e$results$ccg_name),
               date = last_date,
               n_increase = e$results %>%
                  filter(date >= last_date -obs_period + 1 & classification == 'increase') %>%
                  nrow()
            )
         ) %>%
            bind_rows() %>%
            mutate(rank = frankv(n_increase, ties.method='min', order=-1)) %>%
            arrange(rank)
         
         ccg_rank_date <- bind_rows(ccg_rank_date, ccg_rank)
      }
      
      for (sc in selected_ccgs) {
         asmodee_selected_plot[[asmodee_conf]][[sc]] <- list()
         for (i in 1:length(date_range_selected)) {
            last_date <- date_range_selected[i]
            asmodee_selected_plot[[asmodee_conf]][[sc]][[i]] <- plot(
               res_ccg_calls_list[[asmodee_conf]][[i]][[sc]],
               'date', point_size=1, guide=F) +
               ggtitle(format(last_date, '%d %b %Y')) +
               theme(text = element_text(size = 12),
                     axis.text.x = element_text(angle = 45, hjust = 1)) +
               scale_x_date(date_labels = format("%d %b")) +
               geom_vline(xintercept = last_date - obs_period + 1 - 0.6, size=0.75, color='grey80')
         }
      }
      
      ccg_rank_date <- ccg_rank_date %>%
         complete(ccg_name, date, fill=list(n_increase=NA, rank=NA)) %>%
         arrange(ccg_name, date)
      ccg_rank_date_all <- bind_rows(
         ccg_rank_date_all,
         ccg_rank_date %>% mutate(asmodee_conf=asmodee_conf)
      )
      
      
   }
   saveRDS(res_ccg_calls_list, here(data_relative_path, 'res_ccg_calls_list.rds'))
   saveRDS(asmodee_selected_plot, here(data_relative_path, 'asmodee_selected_plot.rds'))
   saveRDS(ccg_rank_date_all, here(data_relative_path, 'ccg_rank_date_all.rds'))
   
    t10 <- Sys.time()
    print('Time elapsed for Leicester:')
    print(t10-t9)
   
} else {
   
   res_ccg_calls_list <- readRDS(here(data_relative_path, 'res_ccg_calls_list.rds'))
   asmodee_selected_plot <- readRDS(here(data_relative_path, 'asmodee_selected_plot.rds'))
   ccg_rank_date_all <- readRDS(here(data_relative_path, 'ccg_rank_date_all.rds'))
   
}

if (plot_ccgs) {
   
   for (asmodee_conf in unique(ccg_rank_date_all$asmodee_conf)) {
      
      ccg_rank_date <- ccg_rank_date_all %>% filter(asmodee_conf==asmodee_conf)
      
      for (sc in selected_ccgs) {
         ccg_name <- CleanCCGName(sc)
         asmodee_selected_plot_combined <- arrangeGrob(
            grobs=asmodee_selected_plot[[asmodee_conf]][[sc]],
            as.table=T,
            ncol=min(length(asmodee_selected_plot[[asmodee_conf]][[sc]]), 7),
            top=textGrob(ccg_name, gp=gpar(fontsize=30))
         )
         ggsave(asmodee_selected_plot_combined,
                filename=here(ccgs_relative_path,
                              paste0('asmodee_selected_', sc, '_', asmodee_conf, '.pdf')),
                width=8*min(length(asmodee_selected_plot[[asmodee_conf]][[sc]]), 7),
                height=6*(1+floor(length(asmodee_selected_plot[[asmodee_conf]][[sc]])/7)), unit='cm')
         ggsave(asmodee_selected_plot_combined,
                filename=here(ccgs_relative_path,
                              paste0('asmodee_selected_', sc, '_', asmodee_conf, '.png')),
                width=8*min(length(asmodee_selected_plot[[asmodee_conf]][[sc]]), 7),
                height=6*(1+floor(length(asmodee_selected_plot[[asmodee_conf]][[sc]])/7)), unit='cm',
                dpi=150)
      }
      
      ccgs_colors <- c(
         vega_standard_palette[1:length(selected_ccgs)],
         rep('grey70', length(unique(ccg_rank_date$ccg_name))-length(selected_ccgs))
      )
      names(ccgs_colors) <- c(
         selected_ccgs,
         unique(ccg_rank_date$ccg_name)[! unique(ccg_rank_date$ccg_name) %in% selected_ccgs]
      )
      ccgs_compare_title <- paste0(
         'How ',
         paste(sapply(selected_ccgs[-length(selected_ccgs)], function (sc)
            paste0('<b style="color:', ccgs_colors[[sc]], '">', CleanCCGName(sc), '</b>')),
            collapse=', '),
         ifelse(length(selected_ccgs)==1, '', ' and '),
         paste0('<b style="color:', ccgs_colors[[tail(selected_ccgs,1)]], '">',
                CleanCCGName(tail(selected_ccgs,1)), '</b>'),
         ifelse(length(selected_ccgs)==1, ' compares', ' compare'),
         ' to other CCG\'s.'
      )
      ccg_selected_nincrease_plot <- ggplot(
         ccg_rank_date %>% filter(!ccg_name %in% selected_ccgs),
         aes(date, n_increase, fill=ccg_name, color=ccg_name)) +
         geom_line(alpha=0.3,
                   position = position_jitter(width=0.05, height=0.1, seed = 1)) +
         geom_point(color='black', shape=21, alpha=0.3,
                    position = position_jitter(width=0.07, height=0.1, seed = 1)) +
         geom_line(data=ccg_rank_date %>% filter(ccg_name %in% selected_ccgs),
                   color='black', size=1.4,
                   position = position_jitter(width=0.07, height=0.1, seed = 2))
      for (ccg in names(lockdown_start)) {
         ccg_selected_nincrease_plot <- ccg_selected_nincrease_plot +
            geom_vline(xintercept=lockdown_start[[ccg]]-0.5, color=ccgs_colors[[ccg]], size=1, alpha=0.75)
      }
      ccg_selected_nincrease_plot <- ccg_selected_nincrease_plot +
         geom_line(data=ccg_rank_date %>% filter(ccg_name %in% selected_ccgs), size=0.9,
                   position = position_jitter(width=0.07, height=0.1, seed = 2)) +
         geom_point(data=ccg_rank_date %>% filter(ccg_name %in% selected_ccgs),
                    color='black', shape=21, size=2,
                    position = position_jitter(width=0.07, height=0.1, seed = 2)) +
         scale_color_manual(values=ccgs_colors) +
         scale_fill_manual(values=ccgs_colors) +
         scale_y_continuous(breaks=0:max(ccg_rank_date$n_increase)) +
         ggtitle(ccgs_compare_title) +
         ylab('number of *increase* days') +
         theme_bw() +
         theme(
            text = element_text(size = 12),
            plot.title = element_markdown(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_markdown(),
            legend.position = 'none',
            panel.grid.minor.y = element_blank()) +
         scale_x_date(date_labels = format("%d %b"))
      
      ggsave(ccg_selected_nincrease_plot,
             filename=here(ccgs_relative_path, paste0('ccg_selected_nincrease_plot_', asmodee_conf, '.pdf')),
             width=20, height=12, unit='cm')
      ggsave(ccg_selected_nincrease_plot,
             filename=here(ccgs_relative_path, paste0('ccg_selected_nincrease_plot_', asmodee_conf, '.png')),
             width=20, height=12, unit='cm', dpi=150)
      
      ccg_selected_rank_plot <- ggplot(
         ccg_rank_date %>% filter(! ccg_name %in% selected_ccgs),
         aes(date, rank, fill=ccg_name, color=ccg_name)) +
         geom_line(alpha=0.3,
                   position = position_jitter(width=0.1, height=1.5, seed = 1)) +
         geom_point(color='black', shape=21, alpha=0.3,
                    position = position_jitter(width=0.1, height=1.5, seed = 1)) +
         geom_line(data=ccg_rank_date %>% filter(ccg_name %in% selected_ccgs),
                   color='black', size=1.4,
                   position = position_jitter(width=0.1, height=1.5, seed = 2))
      for (ccg in names(lockdown_start)) {
         ccg_selected_nincrease_plot <- ccg_selected_nincrease_plot +
            geom_vline(xintercept=lockdown_start[[ccg]]-0.5, color=ccgs_colors[[ccg]])
      }
      ccg_selected_rank_plot <- ccg_selected_rank_plot +
         geom_line(data=ccg_rank_date %>% filter(ccg_name %in% selected_ccgs), size=0.9,
                   position = position_jitter(width=0.1, height=1.5, seed = 2)) +
         geom_point(data=ccg_rank_date %>% filter(ccg_name %in% selected_ccgs),
                    color='black', shape=21, size=2,
                    position = position_jitter(width=0.1, height=1.5, seed = 2)) +
         scale_color_manual(values=ccgs_colors) +
         scale_fill_manual(values=ccgs_colors) +
         scale_y_reverse(breaks = c(1, 30, 60, 90, 120)) +
         theme_bw() +
         theme(text = element_text(size = 12),
               axis.text.x = element_text(angle = 45, hjust = 1),
               legend.position = 'none',
               panel.grid.minor.y = element_blank()) +
         scale_x_date(date_labels = format("%d %b"))
      
      ccg_selected_plot_combined <- arrangeGrob(
         grobs=list(ccg_selected_nincrease_plot + xlab(NULL), ccg_selected_rank_plot),
         as.table=T,
         ncol=1
      )
      ggsave(ccg_selected_plot_combined,
             filename=here(ccgs_relative_path, paste0('ccg_selected_plot_combined_', asmodee_conf, '.pdf')),
             width=25, height=20, unit='cm')
      ggsave(ccg_selected_plot_combined,
             filename=here(ccgs_relative_path, paste0('ccg_selected_plot_combined_', asmodee_conf, '.png')),
             width=25, height=20, unit='cm', dpi=150)
      
      ccg_selected_nincrease_heatmap <- ggplot(ccg_rank_date,
                                               aes(date, CleanCCGName(ccg_name), fill=n_increase)) +
         geom_tile() +
         scale_fill_gradient(low=vega_standard_palette[1], high=vega_standard_palette[2]) +
         labs(y='CCG', fill='number of\nincrease days') +
         theme_bw() +
         theme(panel.grid = element_blank())
      ggsave(ccg_selected_nincrease_heatmap,
             filename=here(ccgs_relative_path,
                           paste0('ccg_selected_nincrease_heatmap_', asmodee_conf, '.pdf')),
             width=25, height=50, unit='cm')
      ggsave(ccg_selected_nincrease_heatmap,
             filename=here(ccgs_relative_path,
                           paste0('ccg_selected_nincrease_heatmap_', asmodee_conf, '.png')),
             width=25, height=50, unit='cm', dpi=150)
      
      ccg_selected_rank_heatmap <- ggplot(ccg_rank_date,
                                          aes(date, CleanCCGName(ccg_name), fill=rank)) +
         geom_tile() +
         scale_fill_gradient(high=vega_standard_palette[1], low=vega_standard_palette[2]) +
         labs(y='CCG', fill='rank') +
         theme_bw() +
         theme(panel.grid = element_blank())
      ggsave(ccg_selected_rank_heatmap,
             filename=here(ccgs_relative_path,
                           paste0('ccg_selected_rank_heatmap_', asmodee_conf, '.pdf')),
             width=25, height=50, unit='cm')
      ggsave(ccg_selected_rank_heatmap,
             filename=here(ccgs_relative_path,
                           paste0('ccg_selected_rank_heatmap_', asmodee_conf, '.png')),
             width=25, height=50, unit='cm', dpi=150)
   }
   
}
ccg_rank_date_all

