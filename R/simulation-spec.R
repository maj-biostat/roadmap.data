library(data.table)

pkg.sim <- new.env()

get_sim_spec <- function(){

  a_s_u <- qlogis(
    array(
      c(0.65, 0.55, 0.60,
        0.75, 0.60, 0.65),
      dim = c(3,2),
      dimnames = list(
        c("early", "late", "chronic"),
        c("knee", "hip")
      )
    )
  )

  # domain membership adjustment
  gamma <- log(
    array(
      c(1,1,1,1,1,0.9,
        1,1,1,1,0.9,0.95,
        1,1,1,1,0.9,0.95),
      dim = c(3,2,3),
      dimnames = list(
        c("a", "b", "c"),
        c("Y", "N"),
        c("early", "late", "chronic")
      )
    )
  )

  # e.g.
  # gamma["c","N","late"]
  # or
  # m <- cbind(c("a", "b", "a", "c"), c("N", "N", "N", "Y"), c("early", "late", "late", "chronic"))
  # gamma[m]


  # effective trt arms

  # surgery

  b_a_late <- c(dair = 0, rev = log(1.5))
  b_a_chronic <- c(one = 0, two = log(1.8))

  # duration
  b_b1_late_one <- c(w06p1 = 0, w12p1 = log(2))
  b_b2_late_two <- c(d07p2 = 0, w12p2 = log(1.4))
  b_b1_chronic_one <- c(w06p1 = 0, w12p1 = log(1.6))
  b_b2_chronic_two <- c(d07p2 = 0, w12p2 = log(1.7))

  # type
  b_c <- c(norif = 0, rif = log(1.3))

  pkg.sim$spec <- list(
    a_s_u = a_s_u,
    gamma = gamma,
    b_a_late = b_a_late,
    b_a_chronic = b_a_chronic,
    b_b1_late_one = b_b1_late_one,
    b_b2_late_two = b_b2_late_two,
    b_b1_chronic_one = b_b1_chronic_one,
    b_b2_chronic_two = b_b2_chronic_two,
    b_c = b_c
  )

  pkg.sim$spec
}
