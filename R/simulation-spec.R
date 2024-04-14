library(data.table)

pkg.sim <- new.env()

get_sim_spec <- function(){


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
    "erx" = -0.1, "erx-r1" = -0.05, "erx-r2" = 0.05,
    "r1" = 0.2, "r2" = 0.09,
    "r1d" = 0.3, "r2d" = 0.1,
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


