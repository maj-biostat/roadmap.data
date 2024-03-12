library(data.table)
library(poisson)

N = 1000
pop_spec = NULL
sim_spec = NULL
idx_s = 1



# Design contains variables:
# l silo membership (early, late, chronic)
# l1 indicator for late silo membership
# l2 indicator for late silo membership
# j indicator for hip joint
# er indicator for reveal of surgical domain
# ed indicator for reveal of duration domain
# ef indicator for reveal of choice domain
# erx indicator for non-reveal of surgical domain
# edx efx
# r random assignment dair/rev
# sr pref dair/one/two for unit i
# sra pref one/two for unit i
# ra allocation of surgical procedure (dair, one, two)
# ic indicator of cross over
# rp surgery received dair/rev
# srp selection of actual procedure
# srp2 indicator of two stage procedure performed
# d duration random assignment
# f ab choice random assignment

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

  j <- rbinom(N, 1, pop_spec$r_joint[l+1, 2])
  # reveal for late only
  # minic the siutation where a small number who never get revealed even if they
  # had late-stage infection.
  # you can leave this out, but the model spec will need to be updated
  # because there will be linearly depenedent cols in the design matrix, e.g.
  # l1 == er

  # if rand has been shut off then do not enrol any more into surgical domain
  if(all(is.na(unlist(pop_spec$r_a)))){
    er <- rep(0, N)
  } else{
    # only late silo can enter surgical domain
    er <- as.numeric(l == 1)
    # 2% are never revealed (for whatever reason)
    i_rec <- as.logical(rbinom(er[l==1], 1, 0.02))
    er[l==1][i_rec] <- 0
  }


  # randomise surgery - silo specific allocation probs.
  r <- rep(NA, N)
  for(i in 1:length(pop_spec$r_a)){
    z <- pop_spec$r_a[[i]]
    if(all(!is.na(z))){
      r[l==(i-1)] <- sample(0:(length(z)-1), sum(l==(i-1)), TRUE, z)
    } else {
      r[l==(i-1)] <- 0
    }
  }
  # fix those that were not randomised to control. in practice r is undef
  # for these units
  r[er == 0] <- 0

  # assume clinician states a preferred surgery procedure for each pt
  # (approx) 70% chance of clinician choosing two-stage if pt is rand to revision
  sr <- rep(NA, N)
  for(i in 1:length(pop_spec$r_a_q)){
    z <- pop_spec$r_a_q[[i]]
    if(all(!is.na(z))){
      # preference out of all of dair, one, two (for those rand to dair)
      r_0 <- (r == 0) & (l == i - 1)
      sr[r_0] <- sample(0:(length(z)-1), sum(r_0), TRUE, z)
      # preference out of one, two (for those that were rand to rev)
      r_1 <- (r == 1) & (l == i - 1)
      sr[r_1] <- sample(1:(length(z)-1), sum(r_1), TRUE, z[2:length(z)]/sum(z[2:length(z)]))
    } else {
      sr[(l == i - 1)] <- 0
    }
  }

  # pref towards two-stage, assuming revision
  sra <- as.numeric(sr == 2)
  # determine allocation of surgery type
  ra <- (1-er) * sr + er * r * (sr)

  # 0.5% of the allocated treatments switch to a different surg type (constrain to min of 1 unit)
  srp <- ra
  ic <- rbinom(N, 1, 0.005)
  if(all(ic == 0)){ ic[sample(length(ic), 1)] <- 1 } # add at least one cross over (there is always one)
  # allows for crossover to any surgical type - these probabilities should really be silo specific
  srp[as.logical(ic)] <- sample(0:2, sum(ic), replace = T, prob = c(0.1, 0.3, 0.6))
  # was the procedure that was performed dair or revision?
  rp <- as.numeric(srp %in% 1:2)

  # reveal for duration

  # reveal only occurs for the duration that the effects are being evaluated
  # if the question of duration is answered in the cohort recv one-stage
  # (indicated by all randomisation probs being NA for that cohort) then
  # they no longer contribute to the randomised comparison for the duration
  # domain.

  dtmp <- data.table(rp, srp, ed = rep(0, N))
  if(all(!is.na(pop_spec$r_b$one))){
    dtmp[rp == 1 & srp == 1, ed := 1]
  }
  if(all(!is.na(pop_spec$r_b$two))){
    dtmp[rp == 1 & srp == 2, ed := 1]
  }
  ed <- dtmp$ed

  # rand to long (0), short (1) for one-stage
  # rand to short (0), long (1) for two-stage

  d <- rep(NA, N)
  for(i in 1:length(pop_spec$r_b)){
    z <- pop_spec$r_b[[i]]
    if(all(!is.na(z))){
      d[srp==(i-1)] <- sample(0:(length(z)-1), sum(srp==(i-1)), TRUE, z)
    } else {
      d[srp==(i-1)] <- 0
    }
  }

  # reveal for ab choice

  # again, reveal only occurs for the duration that the effects are being
  # evaluated if the question of choice is decided
  # (indicated by all randomisation probs being NA for rif) then
  # they no further enrolments happen for ab choice.

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
    r, sr, sra, ra,
    ic, rp,
    srp,
    srp1 = as.numeric(srp == 1),
    srp2 = as.numeric(srp == 2),
    d,
    f
  )

  D[, id := idx_s:(N+idx_s - 1)]
  setcolorder(D, "id")

  D
}

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
        b["erx"] * d$erx +
        (b["r1"] * d$r + b["r2"] * d$r * d$srp2) * d$er +
        b["edx"] * d$edx +
        (b["r1d"] * d$rp * d$d * d$srp1 + b["r2d"] * d$rp * d$d * d$srp2) * d$ed +
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

get_stan_data <- function(d){

  d_s <- d[, .(y = sum(y), n = .N),
           keyby = .(l1, l2, er, ed, ef, r, rp, srp, srp2, d, f)]

  ld <- list(
    N = nrow(d_s), y = d_s$y, n = d_s$n,
    l1 = d_s$l1, l2 = d_s$l2,
    er = d_s$er, ed = d_s$ed, ef = d_s$ef,
    r = d_s$r, d = d_s$d, f = d_s$f,
    rp = d_s$rp,
    srp1 = as.numeric(d_s$srp == 1),
    srp2 = as.numeric(d_s$srp == 2),
    pri_m_sd = rep(1, 2),
    pri_b_sd = rep(1, 8),
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



test_get_design <- function(){


  source("R/simulation-spec.R")
  source("R/population-spec.R")

  # default setting
  pop_spec <- get_pop_spec()
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] > 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp2 == 0, .N] > 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] > 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] > 0)

  # stop rand to surg
  pop_spec <- get_pop_spec()
  pop_spec$r_a$late['dair'] <- NA
  pop_spec$r_a$late['rev'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] == 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp2 == 0, .N] > 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] > 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] > 0)

  # d1[, .N, keyby = .(ed, srp, d)]

  # stop rand to duration (one-stage)
  pop_spec <- get_pop_spec()
  pop_spec$r_b$one['long'] <- NA
  pop_spec$r_b$one['short'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] > 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp2 == 0, .N] == 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] > 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] > 0)

  # stop rand to duration (two-stage)
  pop_spec <- get_pop_spec()
  pop_spec$r_b$two['long'] <- NA
  pop_spec$r_b$two['short'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] > 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp2 == 0, .N] > 0)
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
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp2 == 0, .N] == 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] == 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] > 0)

  # stop rand to choice
  pop_spec <- get_pop_spec()
  pop_spec$r_c['norif'] <- NA
  pop_spec$r_c['rif'] <- NA
  d1 <- get_design(N = 10000, pop_spec)
  stopifnot("no surgery reveal" = d1[er == 1, .N] > 0)
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp2 == 0, .N] > 0)
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
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp2 == 0, .N] == 0)
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
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp2 == 0, .N] > 0)
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
  stopifnot("no duration reveal (one-stage)" = d1[ed == 1 & srp2 == 0, .N] == 0)
  stopifnot("no duration reveal (two-stage)" = d1[ed == 1 & srp2 == 1, .N] == 0)
  stopifnot("no surgery reveal" = d1[ef == 1, .N] == 0)



}


