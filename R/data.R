library(data.table)
library(poisson)

N = 1000
pop_spec = NULL
sim_spec = NULL
idx_s = 1



# Fields in data output:
# id: pt id
# l: strata 0 early, 1 late, 2 chronic
# l1: strata indicator for late
# l2: strata indicator for chronic
# j: joint indicator knee/hip
# er: revealed indicator surgery
# ed: revealed indicator duration
# ef: revealed indicator choice
# erx: 1 - revealed = non-revealed indicator for surgery
# edx: 1 - revealed = non-revealed indicator for duration
# efx: 1 - revealed = non-revealed indicator for choice
# r: randomisation for surgery domain, dair vs rev hard-coded restriction to late
# sr: preferred surgery at elicited at baseline (0 dair, 1 one-stage, 2 two-stage)
# sra: indicator derived from sr for preference for two-stage
# ra: allocated surgical approach accounting for whether rand in surg or not
# ic: indicator for treatment switch (what was planned was not received).
# rp: indicator of performed surgical approach (0 dair, 1 rev)
# srp: performed surgical approach (0 dair, 1 one-stage, 2 two-stage)
# srp1: indicator for one-stage performed
# srp2: indicator for two-stage performed
# Duration domain depends on what surgery was received, NOT what was planned.
# Because of the questions of interest:
# For one-stage long (0), short (1)
# For two-stage short (0), long (1)
# d: indicator for short/long duration
# f: indicator for choice (0 no-rif, 1 rif)
# t0: enrolment time
# eta_y: log-odds treatment success
# p_y: pr treatment success
# y: observed outcome (0 fail, 1 success)

# N = 1e6
# idx_s = 1
# set.seed(1)
# pop_spec <- get_pop_spec()
get_design <- function(N = 2500, pop_spec = NULL, idx_s = 1){

  if(is.null(pop_spec)){
    pop_spec <- get_pop_spec()
  }

  # silo (l)
  l <- sample(0:2, N, replace = T, prob = pop_spec$r_silo )
  # convenience indicators for late and chronic
  l1 = as.numeric(l == 1)
  l2 = as.numeric(l == 2)

  # joint is matrix, rows are silo, cols are knee, hip. j = 1 indicates hip
  j <- rbinom(N, 1, pop_spec$r_joint[l+1, 2])
  # reveal for late only with negligible numbers not entering in for randomised
  # treatment

  # if rand has been shut off then do not enrol any more into surgical domain
  if(all(is.na(unlist(pop_spec$r_a)))){
    # no surgical allocation revealed.
    er <- rep(0, N)
  } else{
    # only late silo can enter surgical domain
    er <- as.numeric(l == 1)
    # but say 2% are never revealed in this silo (for whatever reason)
    i_rec <- as.logical(rbinom(er[l==1], 1, 0.02))
    er[l==1][i_rec] <- 0
  }

  # randomisation assignment for surgery - silo specific allocation probs.
  r <- rep(NA, N)
  for(i in 1:length(pop_spec$r_a)){
    # r_a has allocation probs for early, late and chronic silo
    # pick the current silo
    z <- pop_spec$r_a[[i]]
    # if the elements within the silo are all non-na then
    if(all(!is.na(z))){
      # allocate the surgical trt based on the the number of trt for the strata
      r[l==(i-1)] <- sample(0:(length(z)-1), sum(l==(i-1)), TRUE, z)
    } else {
      # otherwise just set to zero. er will be set to zero so these will be
      # ignored anyway
      r[l==(i-1)] <- 0
    }
  }
  # fix those that are in the late silo but were not revealed, see above
  r[er == 0] <- 0

  # Just take r_a_q to indicate the surgery type that actually happens
  # as preference is not really required any more in the sims.
  # For each silo r_a_q indicates the prob of each surg type (dair, one, two)
  # We still will record preference for the trial.
  dtmp <- data.table(cbind(l, er, r, srp = rep(NA, N)))
  # For each silo.
  for(i in 1:length(pop_spec$r_a_q)){
    z <- pop_spec$r_a_q[[i]]

    # non-randomised trt just gets whatever r_a_q implies
    dtmp[l == i - 1 & er == 0, srp := sample(0:(length(z)-1), .N, TRUE, z)]
    # randomised trt to dair gets dair
    dtmp[l == i - 1 & er == 1 & r == 0, srp := 0]
    # randomised trt to rev gets one/two based on normalised probs
    dtmp[l == i - 1 & er == 1 & r == 1,
         srp := sample(1:(length(z)-1), .N, TRUE, z[2:length(z)]/sum(z[2:length(z)]))]

  }

  # for kicks, set 2% of the late silo to not revealed for randomised surgery
  # and set them to receive dair
  # ic <- rbinom(dtmp[l == 1, .N], 1, 0.02)
  # if(all(ic == 0)){ ic[1] <- 1}
  # dtmp[l == 1, ic := ic]
  # dtmp[l != 1, ic := 0]
  # dtmp[ic == 1, `:=`(er = 0, r = 0, srp = 0)]

  srp <- dtmp$srp

  # rp is now redundant.....?

  # was the procedure that was performed dair or revision?
  # rp <- as.numeric(srp %in% 1:2)

  # reveal for duration

  # Reveal only occurs for the period over which the effects are being evaluated.
  # Once the quest has been answered we stop revealling the randomisation.
  # Presumably we also stop randomising.
  dtmp <- data.table(srp, ed = rep(0, N))
  if(all(!is.na(pop_spec$r_b$one))){
    dtmp[srp == 1, ed := 1]
  }
  if(all(!is.na(pop_spec$r_b$two))){
    dtmp[srp == 2, ed := 1]
  }

  ed <- dtmp$ed

  # rand to long (0), short (1) for one-stage
  # rand to short (0), long (1) for two-stage

  dtmp <- data.table(cbind(srp, ed, d = rep(NA, N)))
  # d <- rep(NA, N)
  # For each surgery type
  for(i in 1:length(pop_spec$r_b)){
    z <- pop_spec$r_b[[i]]
    if(all(!is.na(z))){
      # d[srp==(i-1)] <- sample(0:(length(z)-1), sum(srp==(i-1)), TRUE, z)
      dtmp[srp == (i - 1), d := sample(0:(length(z)-1), .N, TRUE, z)]
    } else {
      dtmp[srp == (i - 1), d := 0]
    }
  }

  d <- dtmp$d

  # reveal for ab choice

  # Again, reveal only occurs for the period that effects are being
  # evaluated. Once the question is answered we reveal no more pts to rand trt.

  if(all(is.na(pop_spec$r_c))){
    ef <- rep(0, N)
    f <- rep(0, N)
  } else {
    # 60% reveal ab choice
    ef <- rbinom(N, 1, 0.6)
    f <- as.numeric((ef == 1) * rbinom(N, 1, pop_spec$r_c))
  }

  D <- data.table(
    l, l1, l2, j,
    er, ed, ef,
    erx = 1-er, edx = 1-ed, efx=1-ef,
    r,
    # sr, sra, ra,
    # ic,
    # rp,
    srp,
    srp0 = as.numeric(srp == 0),
    srp1 = as.numeric(srp == 1),
    srp2 = as.numeric(srp == 2),
    d,
    f
  )

  D[, id := idx_s:(N+idx_s - 1)]
  setcolorder(D, "id")

  D
}

# Generates trial data for cohort size specified by N along with pop_spec and
# sim_spec (if provided otherwise defaults used) according to linear predictor
# specified in g.
# Fields in data output:
# id: pt id
# l: strata 0 early, 1 late, 2 chronic
# l1: strata indicator for late
# l2: strata indicator for chronic
# j: joint indicator knee/hip
# er: revealed indicator surgery
# ed: revealed indicator duration
# ef: revealed indicator choice
# erx: 1 - revealed = non-revealed indicator for surgery
# edx: 1 - revealed = non-revealed indicator for duration
# efx: 1 - revealed = non-revealed indicator for choice
# r: randomisation for surgery domain, dair vs rev hard-coded restriction to late
# sr: preferred surgery at elicited at baseline (0 dair, 1 one-stage, 2 two-stage)
# sra: indicator derived from sr for preference for two-stage
# ra: allocated surgical approach accounting for whether rand in surg or not
# ic: indicator for treatment switch (what was planned was not received).
# rp: indicator of performed surgical approach (0 dair, 1 rev)
# srp: performed surgical approach (0 dair, 1 one-stage, 2 two-stage)
# srp1: indicator for one-stage performed
# srp2: indicator for two-stage performed
# Duration domain depends on what surgery was received, NOT what was planned.
# Because of the questions of interest:
# For one-stage long (0), short (1)
# For two-stage short (0), long (1)
# d: indicator for short/long duration
# f: indicator for choice (0 no-rif, 1 rif)
# t0: enrolment time
# eta_y: log-odds treatment success
# p_y: pr treatment success
# y: observed outcome (0 fail, 1 success)
get_trial_data <- function(
    N = 100000,
    pop_spec = NULL,
    sim_spec = NULL,
    idx_s = 1,
    entry_times = T,
    g = function(d, sim_spec){

      a0 <- sim_spec$a0
      m <- sim_spec$m
      b <- sim_spec$b

      eta <- a0 +
        m["l1"] * d$l1 + m["l2"] * d$l2 +
        (b["erx"] + b["erx-r1"] * d$srp1 + b["erx-r2"] * d$srp2) * d$erx +
        # move to separation of effects for duration based on
        # one-stage and two-stage rather than a linear combination of
        # terms as was done earlier on in develoment.
        # this will necessitate an update to the decision rules
        # and possibly other elements such as simulation par settings
        (b["r1"] * d$srp1 + b["r2"] * d$srp2) * d$r * d$er +

        # exclude this for now as leads to an estimation problem due to
        # collinearity in the data
        # b["edx"] * d$edx +

        # move to separation of effects for duration based on
        # one-stage and two-stage rather than a linear combination of
        # terms as was done earlier on in develoment.

        # 2024-04-03 - think that the d$rp is not required because if
        # srp1 = 1 or srp2 = 1 then revision occurred otherwise dair occurred
        # so you do not need rp.
        # (b["r1d"] * d$rp * d$d * d$srp1 + b["r2d"] * d$rp * d$d * d$srp2) * d$ed +

        (b["r1d"] * d$d * d$srp1 + b["r2d"] * d$d * d$srp2) * d$ed +
        b["efx"] * d$efx +
        b["f"] * d$f * d$ef

      eta
    }
){

  if(is.null(pop_spec)){
    pop_spec <- get_pop_spec()
  }

  if(is.null(sim_spec)){
    sim_spec <- get_sim_spec()
  }

  d <- get_design(N, pop_spec, idx_s)

  if(entry_times == T){
    d[, t0 := get_enrol_time(.N)]
  } else {
    d[, t0 := rep(0.0, .N)]
  }

  d[, eta_y := g(.SD, sim_spec)]
  d[, p_y := plogis(eta_y)]
  d[, y := rbinom(.N, 1, p_y)]

  setkey(d, id)

  list(
    d = d
  )
}

# Wraps up data (expected to have the structure as generated from
# get_trial_data) and puts it into list suitable for stan model.
# Clearly, explicit dependency on model spec so probably doesn't belong here
# and should be refactored out to main sim code.
get_stan_data <- function(d){

  d_s <- d[, .(y = sum(y), n = .N),
           keyby = .(l1, l2, er, r, srp, srp0, srp1, srp2, ed, d, ef, f)]

  ld <- list(
    N = nrow(d_s), y = d_s$y, n = d_s$n,
    l1 = d_s$l1, l2 = d_s$l2,
    er = d_s$er, ed = d_s$ed, ef = d_s$ef,
    r = d_s$r, d = d_s$d, f = d_s$f,
    srp0 = as.numeric(d_s$srp == 0),
    srp1 = as.numeric(d_s$srp == 1),
    srp2 = as.numeric(d_s$srp == 2),
    pri_m_sd = rep(1, 2),
    pri_b_sd = rep(1, 9),
    prior_only = 0
  )

  list(
    d_s = d_s,
    ld = ld
  )

}

get_enrol_time <- function(N = 2500, lambda = 1.52,
                           rho = function(t) pmin(t/360, 1)){

  c(0, nhpp.event.times(lambda, N - 1, rho))
}

get_design_opts <- function(){

  d_x <- CJ(
    l = c("e", "l", "c"),
    er = c("rand", "nonrand"),
    r = c("dair", "rev"),
    sr = c("dair", "one", "two"),
    ed = c("rand", "nonrand"),
    d = c("short", "long"),
    ef = c("rand", "nonrand"),
    f = c("norif", "rif")
  )

  d_x <- d_x[!(er == "nonrand" & ed == "nonrand" & ef == "nonrand")]
  d_x <- d_x[!(l == "e" & er == "rand")]
  d_x <- d_x[!(l == "c" & er == "rand")]
  d_x <- d_x[!(l == "l" & er == "nonrand")]

  d_x <- d_x[!(l == "e" & sr == "two")]

  d_x <- d_x[!(r == "dair" & sr == "one")]
  d_x <- d_x[!(r == "dair" & sr == "two")]
  d_x <- d_x[!(r == "rev" & sr == "dair")]
  d_x <- d_x[!(r == "dair" & ed == "rand")]

  d_x <- d_x[!(ed == "nonrand" & r == "rev")]
  d_x <- d_x[!(ef == "nonrand" & f == "rif")]

  fwrite(d_x, "x.csv")



}

test_get_design <- function(){


  source("R/simulation-spec.R")
  source("R/population-spec.R")

  # default setting
  pop_spec <- get_pop_spec()
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] > 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp1 == 1, .N] > 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] > 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] > 0)

  # stop rand to surg
  pop_spec <- get_pop_spec()
  pop_spec$r_a$late['dair'] <- NA
  pop_spec$r_a$late['rev'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] == 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp1 == 1, .N] > 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] > 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] > 0)

  # d1[, .N, keyby = .(ed, srp, d)]

  # stop rand to duration (one-stage)
  pop_spec <- get_pop_spec()
  pop_spec$r_b$one['long'] <- NA
  pop_spec$r_b$one['short'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] > 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp1 == 1, .N] == 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] > 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] > 0)

  # stop rand to duration (two-stage)
  pop_spec <- get_pop_spec()
  pop_spec$r_b$two['long'] <- NA
  pop_spec$r_b$two['short'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] > 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp1 == 1, .N] > 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] == 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] > 0)

  # stop rand to duration (one-stage and two-stage)
  pop_spec <- get_pop_spec()
  pop_spec$r_b$one['long'] <- NA
  pop_spec$r_b$one['short'] <- NA
  pop_spec$r_b$two['long'] <- NA
  pop_spec$r_b$two['short'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] > 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp1 == 1, .N] == 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] == 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] > 0)

  # stop rand to choice
  pop_spec <- get_pop_spec()
  pop_spec$r_c['norif'] <- NA
  pop_spec$r_c['rif'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] > 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp1 == 1, .N] > 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] > 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] == 0)

  # stop surg and duration
  pop_spec <- get_pop_spec()
  pop_spec$r_a$late['dair'] <- NA
  pop_spec$r_a$late['rev'] <- NA
  pop_spec$r_b$one['long'] <- NA
  pop_spec$r_b$one['short'] <- NA
  pop_spec$r_b$two['long'] <- NA
  pop_spec$r_b$two['short'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] == 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp1 == 1, .N] == 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] == 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] > 0)

  # stop surg and choice
  pop_spec <- get_pop_spec()
  pop_spec$r_a$late['dair'] <- NA
  pop_spec$r_a$late['rev'] <- NA
  pop_spec$r_c['norif'] <- NA
  pop_spec$r_c['rif'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] == 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp1 == 1, .N] > 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] > 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] == 0)

  # stop duration and choice
  pop_spec <- get_pop_spec()
  pop_spec$r_b$one['long'] <- NA
  pop_spec$r_b$one['short'] <- NA
  pop_spec$r_b$two['long'] <- NA
  pop_spec$r_b$two['short'] <- NA
  pop_spec$r_c['norif'] <- NA
  pop_spec$r_c['rif'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] > 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp1 == 1, .N] == 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] == 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] == 0)



}


