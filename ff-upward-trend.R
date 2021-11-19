library(surveillance)
library(here)

# The Farrington algorithm ignores trend if the estimated case count mean is above past observations

ff_control <- list(
  alpha=0.05,
  b=3,
  w=3,
  limit54=c(0,1),
  noPeriods=7,
  pastWeeksNotIncluded=7,
  glmWarnings=FALSE,
  thresholdMethod='nbPlugin',
  pThresholdTrend=1
)

observ <- c(108, 121, 140, 145, 143, 152, 141, 169, 162, 188, 203, 241, 221, 261, 280, 255, 298,
  349, 356, 421, 352, 455, 458, 482, 499, 547, 576, 565, 611, 648, 700, 726, 789, 824, 894, 906,
  870, 819, 829, 778, 765, 768)

sts <- sts(observ, frequency=7)
ff_out <- farringtonFlexible(sts, control=ff_control)

plot(observ)
plot(ff_out)
print(ff_out@control$trend)

# A new version of the `farringtonFlexible` function written by Michael Höhle:
# - offers the option of removing the safeguard against extrapolating trend through
#   `safeguardAgainstExtrapolation=F`
# - defines bounds on one-side p-values of the prediction interval
source(here('farrington', 'farringtonFlexible.R'))
ff_out2 <- farringtonFlexible(sts, control = append(list(safeguardAgainstExtrapolation=F), ff_control))

plot(observ)
plot(ff_out2)
print(ff_out2@control$trend)

# Test definition of upper bound
ff_out2_ub <- ff_out2@upperbound
ff_out2_ub_pi <- qnbinom(
  1-ff_control$alpha,
  size=ff_out2@control$mu0Vector/(ff_out2@control$phiVector-1),
  mu=ff_out2@control$mu0Vector
)
print(sum(ff_out2_ub==ff_out2_ub_pi)/length(ff_out2_ub))
