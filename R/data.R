library(data.table)
library(poisson)

N = 100000
pop_spec = NULL
sim_spec = NULL

get_design <- function(N = 100000, pop_spec = NULL){

  if(is.null(pop_spec)){
    pop_spec <- get_pop_spec()
  }

  d <- data.table()
  # pt id
  d[, id := 1:N]
  d[, silo := sample(pop_spec$r_silo$silo, size = N, replace = T, prob = pop_spec$r_silo$p)]
  setkey(d, silo)
  setkey(pop_spec$r_joint, silo)
  setkey(pop_spec$r_ea, silo)
  for(z in pop_spec$r_silo$silo){
    d[z, joint :=
        sample(pop_spec$r_joint[z, joint], size = .N, replace = T, prob = pop_spec$r_joint[z, p])
    ]

    # might as well set ea (eligibility for domain a)
    d[z, ea := pop_spec$r_ea[z, rand]]
  }

  # Surgery domain (A)
  setkey(pop_spec$r_a, silo)
  d["early", a := "dair"]
  z <- "late"
  d[z, a := sample(pop_spec$r_a[z, a], size = .N, replace = T, prob = pop_spec$r_a[z, p])]
  z <- "chronic"
  d[z, a := sample(pop_spec$r_a[z, a], size = .N, replace = T, prob = pop_spec$r_a[z, p])]

  # introduce the intended surgical approach (for late stage allocated to revision)
  setkey(d, silo, a)
  setkey(pop_spec$r_a_q, silo, a)
  d[.("late", "rev"), qa :=
      sample(pop_spec$r_a_q[.("late", "rev"), qa], size = .N, replace = T, prob = pop_spec$r_a_q[.("late", "rev"), p])]
  d[is.na(qa), qa := copy(a)]

  # Duration domain (B)
  # eligibility based on a
  d <- merge(d, pop_spec$r_eb[, .(qa, eb = rand)], by = c("qa"), all.x = T)
  setcolorder(d, c("id", "silo", "joint", "ea", "a", "qa","eb"))
  setkey(d, id)

  setkey(d, qa)
  setkey(pop_spec$r_b, qa)
  # all dair get 12wk
  d["dair", b := pop_spec$r_b["dair", unique(b)]]

  setkey(d, qa)
  setkey(pop_spec$r_b, qa)

  z <- "one"
  d[z, b := sample(pop_spec$r_b[z, b], size = .N, replace = T, prob = pop_spec$r_b[z, p])]
  z <- "two"
  d[z, b := sample(pop_spec$r_b[z, b], size = .N, replace = T, prob = pop_spec$r_b[z, p])]
  setkey(d, id)

  # Antibiotic type domain (C)

  # Here I make eligibility random to reflect that some of our cohort will not
  # enter into this domain, irrespective of their randomisation in other domains.
  # Each silo is allowed to have a different proportion of pts entering into
  # domain C.
  setkey(d, silo)
  setkey(pop_spec$r_ec, silo)

  for(z in pop_spec$r_ec$silo){
    d[z, ec :=
        sample(c("Y","N"), size = .N, replace = T, prob = c(pop_spec$r_ec[z, p], 1- pop_spec$r_ec[z, p]))
    ]
  }

  d[ec == "Y",
    c := sample(pop_spec$r_c[, c], size = .N, replace = T, prob = pop_spec$r_c[, p])]
  d[ec == "N", c := "other"]

  setkey(d, id)

  d
}

get_enrol_time <- function(N = 2500, lambda = 1.52,
                           rho = function(t) pmin(t/360, 1)){

  c(0, nhpp.event.times(lambda, N - 1, rho))
}

get_trial_data <- function(N = 100000, pop_spec = NULL, sim_spec = NULL){

  if(is.null(pop_spec)){
    pop_spec <- get_pop_spec()
  }

  if(is.null(sim_spec)){
    sim_spec <- get_sim_spec()
  }

  d <- get_design(N, pop_spec)
  # unique(d[order(silo, joint, ea, a, qa, eb, b, ec, c), .SD, .SDcols = !c("id")])

  d[, t0 := get_enrol_time(.N)]
  # unique(d[order(silo, joint, ea, a, qa, eb, b, ec, c, t0), .SD, .SDcols = !c("id")])

  # initialise log-odds trt success
  d[, alpha_su := sim_spec$a_s_u[cbind(d$silo,d$joint)]]

  # non-membership - todo gamma_a meaningful????
  # e.g. gamma["c","N","late"]
  d[, g_a := sim_spec$gamma[cbind("a", d$ea, d$silo)]]
  d[, g_b := sim_spec$gamma[cbind("b", d$eb, d$silo)]]
  d[, g_c := sim_spec$gamma[cbind("c", d$ec, d$silo)]]

  # introduce treatment effects by silo

  # early stage infection
  # can set b_c for all as is currently pooled.
  setkey(d, ec)
  d[.("Y"), b_c := sim_spec$b_c[c]]
  # set non-members to zero so that we don't have NA kicking about in data
  d[.("N"), b_c := 0]

  # late stage
  # assumes full cohort rand in surg - technically, you do not need gamma_a

  setkey(d, silo, ea)
  d[.("late", "Y"), b_a_late := sim_spec$b_a_late[a]]
  # set the missing to zero so there are no NAs
  d[.("early"), b_a_late := 0]
  d[.("late", "N"), b_a_late := 0]
  d[.("chronic"), b_a_late := 0]

  setkey(d, silo, eb, qa)
  d[.("late", "Y", "one"), b_b1_late_one := sim_spec$b_b1_late_one[b]]
  d[.("late", "Y", "two"), b_b2_late_two := sim_spec$b_b2_late_two[b]]
  # missing values fill in at the end.
  # b_c already set

  # chronic
  # assumes full cohort rand in surg - ie you do not need gamma_a

  setkey(d, silo, ea)
  d[.("chronic", "Y"), b_a_chronic := sim_spec$b_a_chronic[a]]
  # set the missing to zero so there are no NAs
  d[.("early"), b_a_chronic := 0]
  d[.("late"), b_a_chronic := 0]
  d[.("chronic", "N"), b_a_chronic := 0]

  setkey(d, silo, eb, qa)
  d[.("chronic", "Y", "one"), b_b1_chronic_one := sim_spec$b_b1_chronic_one[b]]
  d[.("chronic", "Y", "two"), b_b2_chronic_two := sim_spec$b_b2_chronic_two[b]]
  # b_c was set earlier

  # fill missing (redundant) to zero

  # obvious first - any that are non-members of b (i.e early silo and those rand to dair)
  # need to have all the b params set to zero
  setkey(d, eb)
  d[.("N"), `:=`(b_b1_late_one = 0, b_b2_late_two = 0, b_b1_chronic_one = 0, b_b2_chronic_two = 0)]

  # less obvious
  setkey(d, silo, eb, qa)
  # if late stage revision with two-stage plan then set b domain one stage par to zero
  d[.("late", "Y", "two"), `:=`(b_b1_late_one = 0)]
  # analogous but for one-stage plan
  d[.("late", "Y", "one"), `:=`(b_b2_late_two = 0)]
  # all those in chronic do not contribute to b domain late stage pars (params across silos are indep)
  d[.("chronic", "Y"), `:=`(b_b1_late_one = 0, b_b2_late_two = 0)]

  # same idea as above but for chronic pts
  d[.("chronic", "Y", "two"), `:=`(b_b1_chronic_one = 0)]
  d[.("chronic", "Y", "one"), `:=`(b_b2_chronic_two = 0)]
  d[.("late", "Y"), `:=`(b_b1_chronic_one = 0, b_b2_chronic_two = 0)]

  stopifnot(nrow(d[is.na(b_b1_late_one), .(id, silo, qa, eb, b_b1_late_one, b_b2_late_two, b_b1_chronic_one, b_b2_chronic_two)]) == 0)
  stopifnot(nrow(d[is.na(b_b2_late_two), .(id, silo, qa, eb, b_b1_late_one, b_b2_late_two, b_b1_chronic_one, b_b2_chronic_two)]) == 0)
  stopifnot(nrow(d[is.na(b_b1_chronic_one), .(id, silo, qa, eb, b_b1_late_one, b_b2_late_two, b_b1_chronic_one, b_b2_chronic_two)]) == 0)
  stopifnot(nrow(d[is.na(b_b2_chronic_two), .(id, silo, qa, eb, b_b1_late_one, b_b2_late_two, b_b1_chronic_one, b_b2_chronic_two)]) == 0)

  stopifnot(sum(complete.cases(d)) == nrow(d))

  d[, eta := alpha_su +
      # g_a is redunant (and fixed at zero but I include anyway, more as a
      # reminder than anything else that I may need a g_a in the future)
      g_a + g_b + g_c +
      b_a_late +
      b_a_chronic +
      b_b1_late_one + b_b2_late_two +
      b_b1_chronic_one + b_b2_chronic_two +
      b_c
    ]

  d[, y := rbinom(.N, 1, plogis(eta))]

  setkey(d, id)

  d_i <- get_indexes(d)

  list(
    d = d, # original data
    d_i = d_i # convert to indexes
  )
}


get_indexes <- function(d){

  # replace text with indexes
  d_i <- copy(d[, .SD, .SDcols = !c(
    "alpha_su", "g_a", "g_b", "g_c", "b_c",
    "b_a_late", "b_b1_late_one", "b_b2_late_two",
    "b_a_chronic", "b_b1_chronic_one", "b_b2_chronic_two")])
  setkey(d_i, id)

  # 999 are redundant indexes - never get referenced
  # stan doesn't like NA

  i_early_a <- c(dair = 1)
  i_late_a <- c(dair = 1, rev = 2)
  i_chronic_a <- c(one = 1, two = 2)

  i_early_qa <- c(dair = 1)
  i_late_qa <- c(one = 1, two = 2, dair = 3)
  i_chronic_qa <- c(one = 1, two = 2)

  i_early_b <- c(w12 = 1)
  i_late_b <- c(w06p1 = 1, w12p1 = 2, d07p2 = 1, w12p2 = 2, w12 = 3)
  i_chronic_b <- c(w06p1 = 1, w12p1 = 2, d07p2 = 1, w12p2 = 2)

  i_c <- c(norif = 1, rif = 2, other = 3)

  d_i[, silo := factor(silo, levels = c("early", "late", "chronic"))]
  d_i[, joint := factor(joint, levels = c("knee", "hip"))]

  d_i[, ea := as.integer(ea == "Y")]
  d_i[, eb := as.integer(eb == "Y")]
  d_i[, ec := as.integer(ec == "Y")]

  d_i[silo == "early", a := i_early_a[a]]
  d_i[silo == "late", a := i_late_a[a]]
  d_i[silo == "chronic", a := i_chronic_a[a]]
  d_i[, a := as.integer(a)]

  d_i[silo == "early", qa := i_early_qa[qa]]
  d_i[silo == "late", qa := i_late_qa[qa]]
  d_i[silo == "chronic", qa := i_chronic_qa[qa]]
  d_i[, qa := as.integer(qa)]

  d_i[silo == "early", b := i_early_b[b]]
  d_i[silo == "late", b := i_late_b[b]]
  d_i[silo == "chronic", b := i_chronic_b[b]]
  d_i[, b := as.integer(b)]

  d_i[, c := i_c[c]]
  d_i[, c := as.integer(c)]

  # index for intercept
  d_i[silo == "early" & joint == "knee", su := 1]
  d_i[silo == "early" & joint == "hip", su := 2]
  d_i[silo == "late" & joint == "knee", su := 3]
  d_i[silo == "late" & joint == "hip", su := 4]
  d_i[silo == "chronic" & joint == "knee", su := 5]
  d_i[silo == "chronic" & joint == "hip", su := 6]

  d_i
}
