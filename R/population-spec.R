library(data.table)

pkg.pop <- new.env()

get_pop_spec_new <- function(){

  # proportion allocated to each silo
  r_silo <- c("early" = 0.3, "late" = 0.5, "chronic" = 0.2)
  # proportion with each site of infection by silo
  r_joint <- array(
    c(0.4,0.7,0.5,0.6,0.3,0.5), dim = c(3, 2),
    dimnames = list(c("early", "late", "chronic"), c("knee", "hip")))

  # randomisation probabilities for domain a by silo
  r_a <- list(
    early = NA,
    late = c("dair" = 0.5, "rev" = 0.5),
    chronic = NA
  )

  r_a_q <- list(
    early = c("dair" = 0.9, "one" = 0.1, "two" = 0.0),
    late = c("dair" = 0.2, "one" = 0.24, "two" = 0.56),
    chronic = c("dair" = 0.2, "one" = 0.2, "two" = 0.6)
  )

  r_b <- list(
    dair = NA,
    one = c("long" = 0.5, "short" = 0.5),
    two = c("long" = 0.5, "short" = 0.5)
  )

  r_c <- c("norif" = 0.5, "rif" = 0.5)

  pkg.pop$spec_new <- list(
    r_silo = r_silo,
    r_joint = r_joint,
    r_a = r_a,
    r_a_q = r_a_q,
    r_b = r_b,
    r_c = r_c
  )

  pkg.pop$spec_new
}

# get_pop_spec <- function(){
#   # proportion allocated to each silo
#   r_silo <- data.table(
#     silo = c("early", "late", "chronic"),
#     p = c(0.3, 0.5, 0.2)
#   )
#   # proportion with each site of infection by silo
#   r_joint <- CJ(
#     silo = c("early", "late", "chronic"),
#     joint = c("knee", "hip"), sorted=F
#   )
#   r_joint[, p := c(0.4, 0.6, 0.7, 0.3, 0.5, 0.5)]
#
#   r_joint_new <- array(
#     c(0.4,0.7,0.5,0.6,0.3,0.5), dim = c(3, 2),
#     dimnames = list(c("early", "late", "chronic"), c("knee", "hip")))
#
#
#   # eligibility for domain a
#   r_ea <- data.table(
#     silo = c("early", "late", "chronic"),
#     rand = c("N", "Y", "N")
#   )
#   # randomisation probabilities for domain a by silo
#   r_a <- data.table(
#     silo = c("early", "late","late", "chronic"),
#     a = c(NA, "dair", "rev", NA),
#     p = c(NA, 0.5, 0.5, NA)
#   )
#   # these are the intended surgical approaches based on silo and allocation
#   # it is only applicable for late stage infection pts allocated to revision
#   # who will receive either one-stage or two-stage revision based on the
#   # surgeon's plan, i.e. these are not randomised.
#   r_a_q <- data.table(
#     silo = c(rep("early", 3),rep("late", 3),rep("chronic", 3)),
#     a = c(NA, NA, NA, "rev", "rev", "rev",  NA, NA, NA),
#     # intended/planned surgical approach
#     qa = c("dair", "one", "two", "dair", "one", "two", "dair", "one", "two"),
#     p = c(0.9, 0.1, 0.0,
#           0.3, 0.7, 0.0,
#           0.2, 0.2, 0.6)
#   )
#   # eligibility for domain b
#   r_eb <- data.table(
#     qa = c("dair", "rev", "one", "two"),
#     rand = c("N","Y","Y","Y")
#   )
#   r_b <- data.table(
#     silo = c("early", rep("late", 4), rep("chronic", 4)),
#     qa = c("dair",
#            rep(c("one", "one", "two", "two"), len = 8)),
#     b = c("w12",
#           rep(c("w06p1", "w12p1", "d07p2", "w12p2"), len = 8))
#   )
#   r_b[, p := c(1, rep(0.5, 8))]
#   r_ec <- data.table(
#     silo = c("early", "late","chronic"),
#     rand = c("Y","Y","Y"),
#     p = rep(0.6, 3)
#   )
#   r_c <- data.table(
#     c = c("norif", "rif"), p = rep(0.5, 2)
#   )
#
#   pkg.pop$spec <- list(
#     r_silo = r_silo,
#     r_joint = r_joint,
#     r_ea = r_ea,
#     r_a = r_a,
#     r_a_q = r_a_q,
#     r_eb = r_eb,
#     r_b = r_b,
#     r_ec = r_ec,
#     r_c = r_c
#   )
#
#   pkg.pop$spec
# }


# in case you want to override default spec
set_pop_spec <- function(spec) {
  pkg.pop$spec <- spec
}
