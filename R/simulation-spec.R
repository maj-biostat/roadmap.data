library(data.table)

pkg.sim <- new.env()

# somewhat meaningful names for parameters at relevant indexes
# from stan model.
get_par_effects_mapping <- function(){

  c(
    b_a_l_2 = "b_a_late_rev",
    b_a_c_2 = "b_a_chronic_two",
    b_b1_l_2 = "b_b1_late_one_w12p1",
    b_b2_l_2 = "b_b2_late_two_w12p2",
    b_b1_c_2 = "b_b1_chronic_one_w12p1", # ref level is w06p1
    b_b2_c_2 = "b_b2_chronic_two_w12p2", # ref level is d07p2
    b_c_2 = "b_c_rif"
  )


}

get_sim_spec_new <- function(){


  a0 <- qlogis(0.65)

  # trying to simplify life:
  # just incorporate a shift for silo, don't worry about joint site

  # shift for late, knee
  # m1 <- qlogis(0.55) - a0
  # shift for chronic, knee
  # m2 <- qlogis(0.60) - a0
  # shift for hip
  # m3 <- qlogis(0.75) - a0
  # additional shift for late, hip
  # m4 <- qlogis(0.60) - a0 - m1 - m3
  # additional shift for chronic, hip
  # m5 <- qlogis(0.65) - a0 - m2 - m3
  # sanity check
  # X <- model.matrix(~l*j, data = CJ(j = factor(1:2), l = factor(1:3))[, .(l,j)])

  # a_l_j <- qlogis(
  #   array(
  #     c(plogis(X %*% c(a0, m1, m2, m3, m4, m5))),
  #     dim = c(3,2),
  #     dimnames = list(
  #       c("early", "late", "chronic"),
  #       c("knee", "hip")
  #     )
  #   )
  # )
  # p_l_j <- plogis(a_l_j)

  # the estimates for late and chronic are just weighted combination
  # where the weights correspond to the proportion of joints for each silo

  # shift for late
  m1 <- qlogis(0.57) - a0
  # shift for chronic
  m2 <- qlogis(0.63) - a0

  b <- c(
    "erx" = -0.1, "r1" = 0.2, "r2" = 0.09,
    "edx" = -0.15, "r1d" = 0.3, "r2d" = 0.1,
    "efx" = 0.25, "f" = 0.15
    )

  pkg.sim$spec <- list(
    # a_l_j = a_l_j,
    # p_l_j = p_l_j,
    a0 = a0, m = c("l1" = m1, "l2" = m2),
    b = b
  )

  pkg.sim$spec

}

# get_sim_spec <- function(){
#
#   a_s_u <- qlogis(
#     array(
#       c(0.65, 0.55, 0.60,
#         0.75, 0.60, 0.65),
#       dim = c(3,2),
#       dimnames = list(
#         c("early", "late", "chronic"),
#         c("knee", "hip")
#       )
#     )
#   )
#
#   # domain membership adjustment
#   gamma <- log(
#     array(
#       c(1,1,1,1,1,0.9,
#         1,1,1,1,1,0.9,
#         1,1,1,1,1,0.9),
#       dim = c(3,2,3),
#       dimnames = list(
#         c("a", "b", "c"),
#         c("Y", "N"),
#         c("early", "late", "chronic")
#       )
#     )
#   )
#
#   # e.g.
#   # gamma["c","N","late"]
#   # or
#   # m <- cbind(c("a", "b", "a", "c"), c("N", "N", "N", "Y"), c("early", "late", "late", "chronic"))
#   # gamma[m]
#
#
#   # effective trt arms
#
#   # surgery
#
#   b_a_late <- c(dair = 0, rev = log(1.5))
#   b_a_chronic <- c(one = 0, two = log(1.8))
#
#   # duration
#   b_b1_late_one <- c(w06p1 = 0, w12p1 = log(2))
#   b_b2_late_two <- c(d07p2 = 0, w12p2 = log(1.4))
#   b_b1_chronic_one <- c(w06p1 = 0, w12p1 = log(1.6))
#   b_b2_chronic_two <- c(d07p2 = 0, w12p2 = log(1.7))
#
#   # type
#   b_c <- c(norif = 0, rif = log(1.3))
#
#
#   i_a_s_u <- matrix(1:6, byrow = T, ncol = 2)
#   colnames(i_a_s_u) <- c("knee", "hip")
#   rownames(i_a_s_u) <- c("early", "late", "chronic")
#
#
#   i_a_s_u <- CJ(silo = c("early", "late", "chronic"),
#                 joint = c("knee", "hip"))
#   setkey(i_a_s_u, silo, joint)
#   i_a_s_u[.("early", "knee"), su := 1]
#   i_a_s_u[.("early", "hip"), su := 2]
#   i_a_s_u[.("late", "knee"), su := 3]
#   i_a_s_u[.("late", "hip"), su := 4]
#   i_a_s_u[.("chronic", "knee"), su := 5]
#   i_a_s_u[.("chronic", "hip"), su := 6]
#
#   i_early_a <- c(dair = 1)
#   i_late_a <- c(dair = 1, rev = 2)
#   i_chronic_a <- c(one = 1, two = 2)
#
#   # domain a index mappings
#   i_early_qa <- c(dair = 1)
#   i_late_qa <- c(one = 1, two = 2, dair = 3)
#   i_chronic_qa <- c(one = 1, two = 2)
#
#   # domain b index mappings
#   i_early_b <- c(w12 = 1)
#   i_late_b <- c(w06p1 = 1, w12p1 = 2, d07p2 = 1, w12p2 = 2, w12 = 3)
#   i_chronic_b <- c(w06p1 = 1, w12p1 = 2, d07p2 = 1, w12p2 = 2)
#
#   i_c <- c(norif = 1, rif = 2, other = 3)
#
#
#   pkg.sim$spec <- list(
#     a_s_u = a_s_u,
#     gamma = gamma,
#     b_a_late = b_a_late,
#     b_a_chronic = b_a_chronic,
#     b_b1_late_one = b_b1_late_one,
#     b_b2_late_two = b_b2_late_two,
#     b_b1_chronic_one = b_b1_chronic_one,
#     b_b2_chronic_two = b_b2_chronic_two,
#     b_c = b_c,
#     i_a_s_u = i_a_s_u,
#     i_early_a = i_early_a,
#     i_late_a = i_late_a,
#     i_chronic_a = i_chronic_a,
#     i_early_qa = i_early_qa,
#     i_late_qa = i_late_qa,
#     i_chronic_qa = i_chronic_qa,
#     i_early_b = i_early_b,
#     i_late_b = i_late_b,
#     i_chronic_b = i_chronic_b,
#     i_c = i_c
#   )
#
#   pkg.sim$spec
# }



# get_sim_spec_effects <- function(sim_spec = NULL){
#
#   stopifnot(!is.null(sim_spec))
#
#   data.table(
#     b_a_late_rev = sim_spec$b_a_late["rev"],
#     b_a_chronic_two = sim_spec$b_a_chronic["two"],
#     b_b1_late_one_w12p1 = sim_spec$b_b1_late_one["w12p1"],
#     b_b1_late_two_w12p2 = sim_spec$b_b2_late_two["w12p2"],
#     b_b1_chronic_one_w12p1 = sim_spec$b_b1_chronic_one["w12p1"],
#     b_b2_chronic_two_w12p2 = sim_spec$b_b2_chronic_two["w12p2"],
#     b_c_rif = sim_spec$b_c["rif"]
#   )
# }
