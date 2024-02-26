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
N = 1e6
idx_s = 1
set.seed(1)
pop_spec <- get_pop_spec_new()
get_design_new <- function(N = 2500, pop_spec = NULL, idx_s = 1){

  if(is.null(pop_spec)){
    pop_spec <- get_pop_spec_new()
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
  er <- as.numeric(l == 1)
  i_rec <- as.logical(rbinom(er[l==1], 1, 0.04))
  er[l==1][i_rec] <- 0

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

  # max of 1 unit or 0.5% of the allocated treatments switch to a different surg type
  srp <- ra
  ic <- rbinom(N, 1, 0.005)
  if(all(ic == 0)){ ic[sample(length(ic), 1)] <- 1 } # add at least one cross over (there is always one)
  # allows for crossover to any surgical type - these probabilities should really be silo specific
  srp[as.logical(ic)] <- sample(0:2, sum(ic), replace = T, prob = c(0.1, 0.3, 0.6))
  # was the procedure that was performed dair or revision?
  rp <- as.numeric(srp %in% 1:2)

  # reveal for duration
  ed <- as.numeric(rp == 1)
  # rand to long (0), short (1) based on surgery received

  d <- rep(NA, N)
  for(i in 1:length(pop_spec$r_b)){
    z <- pop_spec$r_b[[i]]
    if(all(!is.na(z))){
      d[srp==(i-1)] <- sample(0:(length(z)-1), sum(srp==(i-1)), TRUE, z)
    } else {
      d[srp==(i-1)] <- 0
    }
  }

  # dtmp <- data.table(id = 1:N, l1, er, r, ra, ic, srp, d)
  # dtmp[er == 1 & r == 0 & l1== 1 & srp == 1 ]
  # dtmp[er == 1 & r == 0 & l1== 1 & srp == 2 ]

  # 60% reveal ab choice
  ef <- rbinom(N, 1, 0.6)
  f <- as.numeric((ef == 1) * rbinom(N, 1, pop_spec$r_c))

  D <- data.table(
    l, l1, l2, j,
    er, ed, ef,
    erx = 1-er, edx = 1-ed, efx=1-ef,
    r, sr, sra, ra,
    ic, rp, srp, srp2 = as.numeric(srp == 2), d,
    f
  )

  D[, id := idx_s:(N+idx_s - 1)]
  setcolorder(D, "id")

  D
}

get_trial_data_new <- function(
    N = 100000,
    pop_spec = NULL,
    sim_spec = NULL,
    idx_s = 1,
    entry_times = T,
    g = function(d, sim_spec){

      a0 <- sim_spec$a0
      m <- sim_spec$m
      b <- sim_spec$b

      # eta <- a0 +
      #   m["l1"] * d$l1 + m["l2"] * d$l2 + m["j"] * d$j +
      #   m["l1j"] * d$l1 * d$j + m["l2j"] * d$l2 * d$j +
      #   b["erx"] * d$erx +
      #   (b["r1"] * d$r + b["r2"] * d$r * d$srp2) * d$er +
      #   b["edx"] * d$edx +
      #   (b["r1d"] * d$rp * d$d + b["r2d"] * d$rp * d$d * d$srp2) * d$ed +
      #   b["efx"] * d$efx +
      #   b["f"] * d$f * d$ef

      eta <- a0 +
        m["l1"] * d$l1 + m["l2"] * d$l2 +
        b["erx"] * d$erx +
        (b["r1"] * d$r + b["r2"] * d$r * d$srp2) * d$er +
        b["edx"] * d$edx +
        (b["r1d"] * d$rp * d$d + b["r2d"] * d$rp * d$d * d$srp2) * d$ed +
        b["efx"] * d$efx +
        b["f"] * d$f * d$ef

      eta
      }
    ){

  if(is.null(pop_spec)){
    pop_spec <- get_pop_spec_new()
  }

  if(is.null(sim_spec)){
    sim_spec <- get_sim_spec_new()
  }

  d <- get_design_new(N, pop_spec, idx_s)

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

get_stan_data_new <- function(d){

  d_s <- d[, .(y = sum(y), n = .N),
           keyby = .(l1, l2, er, ed, ef, r, rp, srp2, d, f)]

  ld <- list(
    N = nrow(d_s), y = d_s$y, n = d_s$n,
    l1 = d_s$l1, l2 = d_s$l2,
    er = d_s$er, ed = d_s$ed, ef = d_s$ef,
    r = d_s$r, d = d_s$d, f = d_s$f,
    rp = d_s$rp, srp2 = d_s$srp2,
    pri_m_sd = rep(1, 2),
    pri_b_sd = rep(1, 8),
    prior_only = 0
  )

  list(
    d_s = d_s,
    ld = ld
  )

}


# Original approach - keep in place until switch over to new version so that
# everything else keeps working.


# get_design <- function(N = 100000, pop_spec = NULL, idx_s = 1){
#
#   if(is.null(pop_spec)){
#     pop_spec <- get_pop_spec()
#   }
#
#   d <- data.table()
#   # pt id
#   d[, id := idx_s:(N+idx_s - 1)]
#   d[, silo := sample(pop_spec$r_silo$silo, size = N, replace = T, prob = pop_spec$r_silo$p)]
#   setkey(d, silo)
#   setkey(pop_spec$r_joint, silo)
#   setkey(pop_spec$r_ea, silo)
#   # for each silo, create joint covariate based on the pop_spec dist
#   # also create entry into domain a
#   for(z in pop_spec$r_silo$silo){
#     d[z, joint :=
#         sample(pop_spec$r_joint[z, joint], size = .N, replace = T, prob = pop_spec$r_joint[z, p])
#     ]
#
#     # might as well set ea (eligibility for domain a)
#     d[z, ea := pop_spec$r_ea[z, rand]]
#   }
#
#   # Surgery domain (A)
#   setkey(pop_spec$r_a, silo)
#   z <- "late"
#
#   d["early", a := "dair"]
#   z <- "late"
#   d[z, a := sample(pop_spec$r_a[z, a], size = .N, replace = T, prob = pop_spec$r_a[z, p])]
#   z <- "chronic"
#   d[z, a := sample(pop_spec$r_a_q[z, qa], size = .N, replace = T, prob = pop_spec$r_a_q[z, p])]
#
#   # introduce the intended surgical approach (for late stage allocated to revision)
#   setkey(d, silo, a)
#   setkey(pop_spec$r_a_q, silo, a)
#   d[.("late", "rev"), qa :=
#       sample(pop_spec$r_a_q[.("late", "rev"), qa], size = .N, replace = T, prob = pop_spec$r_a_q[.("late", "rev"), p])]
#   d[is.na(qa), qa := copy(a)]
#
#   # Duration domain (B)
#   # eligibility based on a
#   d <- merge(d, pop_spec$r_eb[, .(qa, eb = rand)], by = c("qa"), all.x = T)
#   setcolorder(d, c("id", "silo", "joint", "ea", "a", "qa","eb"))
#   setkey(d, id)
#
#   setkey(d, qa)
#   setkey(pop_spec$r_b, qa)
#   # all dair get 12wk
#   d["dair", b := pop_spec$r_b["dair", unique(b)]]
#
#   setkey(d, silo, qa)
#   setkey(pop_spec$r_b, silo, qa)
#
#   # do these separately so that triggers can coordinate silo specific adaptations
#   d[.("late", "one"), b := sample(
#     pop_spec$r_b[.("late", "one"), b], size = .N, replace = T,
#     prob = pop_spec$r_b[.("late", "one"), p])]
#   d[.("late", "two"), b := sample(
#     pop_spec$r_b[.("late", "two"), b], size = .N, replace = T,
#     prob = pop_spec$r_b[.("late", "two"), p])]
#
#   d[.("chronic", "one"), b := sample(
#     pop_spec$r_b[.("chronic", "one"), b], size = .N, replace = T,
#     prob = pop_spec$r_b[.("chronic", "one"), p])]
#   d[.("chronic", "two"), b := sample(
#     pop_spec$r_b[.("chronic", "two"), b], size = .N, replace = T,
#     prob = pop_spec$r_b[.("chronic", "two"), p])]
#
#   setkey(d, id)
#
#   # Antibiotic type domain (C)
#
#   # Here I make eligibility random to reflect that some of our cohort will not
#   # enter into this domain, irrespective of their randomisation in other domains.
#   # Each silo is allowed to have a different proportion of pts entering into
#   # domain C.
#   setkey(d, silo)
#   setkey(pop_spec$r_ec, silo)
#
#   for(z in pop_spec$r_ec$silo){
#     d[z, ec :=
#         sample(c("Y","N"), size = .N, replace = T, prob = c(pop_spec$r_ec[z, p], 1- pop_spec$r_ec[z, p]))
#     ]
#   }
#
#   d[ec == "Y",
#     c := sample(pop_spec$r_c[, c], size = .N, replace = T, prob = pop_spec$r_c[, p])]
#   d[ec == "N", c := "other"]
#
#   setkey(d, id)
#
#   d
# }
#
# get_enrol_time <- function(N = 2500, lambda = 1.52,
#                            rho = function(t) pmin(t/360, 1)){
#
#   c(0, nhpp.event.times(lambda, N - 1, rho))
# }
#
# get_trial_data <- function(N = 100000, pop_spec = NULL, sim_spec = NULL,
#                            idx_s = 1, entry_times = T){
#
#   if(is.null(pop_spec)){
#     pop_spec <- get_pop_spec()
#   }
#
#   if(is.null(sim_spec)){
#     sim_spec <- get_sim_spec()
#   }
#
#   d <- get_design(N, pop_spec, idx_s)
#   # unique(d[order(silo, joint, ea, a, qa, eb, b, ec, c), .SD, .SDcols = !c("id")])
#
#   if(entry_times == T){
#     d[, t0 := get_enrol_time(.N)]
#   } else {
#     d[, t0 := rep(0.0, .N)]
#   }
#
#   # unique(d[order(silo, joint, ea, a, qa, eb, b, ec, c, t0), .SD, .SDcols = !c("id")])
#
#   # initialise log-odds trt success
#   d[, alpha_su := sim_spec$a_s_u[cbind(d$silo,d$joint)]]
#
#   # non-membership - todo gamma_a meaningful????
#   # e.g. gamma["c","N","late"]
#   d[, g_a := sim_spec$gamma[cbind("a", d$ea, d$silo)]]
#   d[, g_b := sim_spec$gamma[cbind("b", d$eb, d$silo)]]
#   d[, g_c := sim_spec$gamma[cbind("c", d$ec, d$silo)]]
#
#   # introduce treatment effects by silo
#
#   # early stage infection
#   # can set b_c for all as is currently pooled.
#   setkey(d, ec)
#   d[.("Y"), b_c := sim_spec$b_c[c]]
#   # set non-members to zero so that we don't have NA kicking about in data
#   d[.("N"), b_c := 0]
#
#   # late stage
#   # assumes full cohort rand in surg - technically, you do not need gamma_a
#
#   setkey(d, silo, ea)
#   d[.("late", "Y"), b_a_late := sim_spec$b_a_late[a]]
#   # d[silo == "late", .(unique(b_a_late), .N), keyby = .(ea, a)]
#   # set the missing to zero so there are no NAs
#   d[.("early"), b_a_late := 0]
#   d[.("late", "N"), b_a_late := 0]
#   d[.("chronic"), b_a_late := 0]
#
#   setkey(d, silo, eb, qa)
#   d[.("late", "Y", "one"), b_b1_late_one := sim_spec$b_b1_late_one[b]]
#   d[.("late", "Y", "two"), b_b2_late_two := sim_spec$b_b2_late_two[b]]
#   # missing values fill in at the end.
#   # b_c already set
#
#   # chronic
#   # assumes full cohort rand in surg - ie you do not need gamma_a
#
#   setkey(d, silo, ea)
#   d[.("chronic", "Y"), b_a_chronic := sim_spec$b_a_chronic[a]]
#   # set the missing to zero so there are no NAs
#   d[.("early"), b_a_chronic := 0]
#   d[.("late"), b_a_chronic := 0]
#   d[.("chronic", "N"), b_a_chronic := 0]
#
#   setkey(d, silo, eb, qa)
#   d[.("chronic", "Y", "one"), b_b1_chronic_one := sim_spec$b_b1_chronic_one[b]]
#   d[.("chronic", "Y", "two"), b_b2_chronic_two := sim_spec$b_b2_chronic_two[b]]
#   # b_c was set earlier
#
#   # fill missing (redundant) to zero
#
#   # obvious first - any that are non-members of b (i.e early silo and those rand to dair)
#   # need to have all the b params set to zero
#   setkey(d, eb)
#   d[.("N"), `:=`(b_b1_late_one = 0, b_b2_late_two = 0, b_b1_chronic_one = 0, b_b2_chronic_two = 0)]
#
#   # less obvious
#   setkey(d, silo, eb, qa)
#   # if late stage revision with two-stage plan then set b domain one stage par to zero
#   d[.("late", "Y", "two"), `:=`(b_b1_late_one = 0)]
#   # analogous but for one-stage plan
#   d[.("late", "Y", "one"), `:=`(b_b2_late_two = 0)]
#   # all those in chronic do not contribute to b domain late stage pars (params across silos are indep)
#   d[.("chronic", "Y"), `:=`(b_b1_late_one = 0, b_b2_late_two = 0)]
#
#   # same idea as above but for chronic pts
#   d[.("chronic", "Y", "two"), `:=`(b_b1_chronic_one = 0)]
#   d[.("chronic", "Y", "one"), `:=`(b_b2_chronic_two = 0)]
#   d[.("late", "Y"), `:=`(b_b1_chronic_one = 0, b_b2_chronic_two = 0)]
#
#   stopifnot(nrow(d[is.na(b_b1_late_one), .(id, silo, qa, eb, b_b1_late_one, b_b2_late_two, b_b1_chronic_one, b_b2_chronic_two)]) == 0)
#   stopifnot(nrow(d[is.na(b_b2_late_two), .(id, silo, qa, eb, b_b1_late_one, b_b2_late_two, b_b1_chronic_one, b_b2_chronic_two)]) == 0)
#   stopifnot(nrow(d[is.na(b_b1_chronic_one), .(id, silo, qa, eb, b_b1_late_one, b_b2_late_two, b_b1_chronic_one, b_b2_chronic_two)]) == 0)
#   stopifnot(nrow(d[is.na(b_b2_chronic_two), .(id, silo, qa, eb, b_b1_late_one, b_b2_late_two, b_b1_chronic_one, b_b2_chronic_two)]) == 0)
#
#   stopifnot(sum(complete.cases(d)) == nrow(d))
#
#   # d[silo == "late", .(.N), keyby = .(
#   #   silo, joint, ea, a, qa, eb, b, ec, c,
#   #   eta,
#   #   b_a_late, b_a_chronic,
#   #   b_b1_late_one, b_b2_late_two, b_b1_chronic_one, b_b2_chronic_two,
#   #   b_c)]
#
#
#   d[, eta := alpha_su +
#       # g_a is redunant (and fixed at zero but I include anyway, more as a
#       # reminder than anything else that I may need a g_a in the future)
#       g_a + g_b + g_c +
#       b_a_late +
#       b_a_chronic +
#       b_b1_late_one + b_b2_late_two +
#       b_b1_chronic_one + b_b2_chronic_two +
#       b_c
#     ]
#
#   d[, y := rbinom(.N, 1, plogis(eta))]
#
#   setkey(d, id)
#
#   d_i <- get_indexes(d, sim_spec)
#
#   list(
#     d = d, # original data
#     d_i = d_i # convert to indexes
#   )
# }
#
#
# get_indexes <- function(d, sim_spec = NULL){
#
#   if(is.null(sim_spec)){
#     sim_spec <- get_sim_spec()
#   }
#
#
#   # replace text with indexes
#   d_i <- copy(d[, .SD, .SDcols = !c(
#     "alpha_su", "g_a", "g_b", "g_c", "b_c",
#     "b_a_late", "b_b1_late_one", "b_b2_late_two",
#     "b_a_chronic", "b_b1_chronic_one", "b_b2_chronic_two")])
#   setkey(d_i, id)
#
#   d_i[, silo := factor(silo, levels = c("early", "late", "chronic"))]
#   d_i[, joint := factor(joint, levels = c("knee", "hip"))]
#
#   d_i[, ea := as.integer(ea == "Y")]
#   d_i[, eb := as.integer(eb == "Y")]
#   d_i[, ec := as.integer(ec == "Y")]
#
#   d_i[silo == "early", a := sim_spec$i_early_a[a]]
#   d_i[silo == "late", a := sim_spec$i_late_a[a]]
#   d_i[silo == "chronic", a := sim_spec$i_chronic_a[a]]
#   d_i[, a := as.integer(a)]
#
#   d_i[silo == "early", qa := sim_spec$i_early_qa[qa]]
#   d_i[silo == "late", qa := sim_spec$i_late_qa[qa]]
#   d_i[silo == "chronic", qa := sim_spec$i_chronic_qa[qa]]
#   d_i[, qa := as.integer(qa)]
#
#   d_i[silo == "early", b := sim_spec$i_early_b[b]]
#   d_i[silo == "late", b := sim_spec$i_late_b[b]]
#   d_i[silo == "chronic", b := sim_spec$i_chronic_b[b]]
#   d_i[, b := as.integer(b)]
#
#   d_i[, c := sim_spec$i_c[c]]
#   d_i[, c := as.integer(c)]
#
#   # index for intercept
#
#   d_i <- merge(d_i, sim_spec$i_a_s_u, by = c("silo", "joint"))
#   setkey(d_i, id)
#
#   d_i
# }
#
# get_stan_data <- function(d_i){
#
#   # ensures that d_b is ordered in the same way as ld so that the indexes
#   # for eta in the generated quantities block line up.
#   d_i[, silo := factor(silo, levels = c("early", "late", "chronic"))]
#   d_i[, joint := factor(joint, levels = c("knee", "hip"))]
#
#   d_b <- d_i[, .(y = sum(y), n = .N), keyby = .(silo, joint, su, ea, a, qa, eb, b, ec, c, eta)]
#
#   ld <- list(
#     N_e = d_b[silo == "early", .N],
#     e_su = d_b[silo == "early", su],
#     e_y = d_b[silo == "early", y],
#     e_n = d_b[silo == "early", n],
#     e_ec = d_b[silo == "early", ec],
#     e_ecp = d_b[silo == "early", 1-ec],
#     e_c = d_b[silo == "early", c],
#
#     N_l = d_b[silo == "late", .N],
#     l_su = d_b[silo == "late", su],
#     l_y = d_b[silo == "late", y],
#     l_n = d_b[silo == "late", n],
#     l_ec = d_b[silo == "late", ec],
#     l_ecp = d_b[silo == "late", 1-ec],
#     l_c = d_b[silo == "late", c],
#     l_ea = d_b[silo == "late", ea],
#     l_eap = d_b[silo == "late", 1-ea],
#     l_a = d_b[silo == "late", a],
#     # below a indicates revision and plan indicates one-stage
#     l_eb1 = d_b[silo == "late", as.integer(a == 2 & qa == 1)],
#     # below a indicates revision and plan indicates two-stage
#     l_eb2 = d_b[silo == "late", as.integer(a == 2 & qa == 2)],
#     l_ebp = d_b[silo == "late", 1-eb],
#     l_b = d_b[silo == "late", b],
#
#     # chronic silo
#     N_c = d_b[silo == "chronic", .N],
#     c_su = d_b[silo == "chronic", su],
#     c_y = d_b[silo == "chronic", y],
#     c_n = d_b[silo == "chronic", n],
#     # domain c randomisation/membership
#     c_ec = d_b[silo == "chronic", ec],
#     # domain c non-randomisation/non-membership
#     c_ecp = d_b[silo == "chronic", 1-ec],
#     # domain c allocation
#     c_c = d_b[silo == "chronic", c],
#     # domain a allocation
#     c_ea = d_b[silo == "chronic", ea],
#     c_eap = d_b[silo == "chronic", 1-ea],
#     c_a = d_b[silo == "chronic", a],
#     # domain b randomisation/membership for one-stage pt
#     c_eb1 = d_b[silo == "chronic", as.integer(a == 1)],
#     # domain b randomisation/membership for two-stage pt
#     c_eb2 = d_b[silo == "chronic", as.integer(a == 2)],
#     # domain b non-randomisation/non-membership
#     c_ebp = d_b[silo == "chronic", 1-eb],
#     c_b = d_b[silo == "chronic", b],
#
#     # default are standard normal priors, just set sd:
#     pri_sig_b_c = 1.0,
#     pri_sig_a_l = 1.0,
#     pri_sig_b1_l = 1.0,
#     pri_sig_b2_l = 1.0,
#     pri_sig_a_c = 1.0,
#     pri_sig_b1_c = 1.0,
#     pri_sig_b2_c = 1.0
#
#   )
#
#   list(
#     d_b = d_b,
#     ld = ld
#   )
#
# }
#




