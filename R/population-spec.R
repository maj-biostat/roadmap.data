library(data.table)

pkg.pop <- new.env()


get_pop_spec <- function(){
  r_silo <- data.table(
    silo = c("early", "late", "chronic"),
    p = c(0.3, 0.5, 0.2)
  )
  r_joint <- CJ(
    silo = c("early", "late", "chronic"),
    joint = c("knee", "hip"), sorted=F
  )
  r_joint[, p := c(0.4, 0.6, 0.3, 0.7, 0.5, 0.5)]
  # eligibility for domain a
  r_ea <- data.table(
    silo = c("early", "late", "chronic"),
    rand = c("N", "Y", "Y")
  )
  r_a <- data.table(
    silo = c("early", "late","late", "chronic", "chronic"),
    a = c("dair", "dair", "rev", "one", "two"),
    p = c(NA, 0.5, 0.5, 0.5, 0.5)
  )
  # these are the intended surgical approaches based on silo and allocation
  # it is only applicable for late stage infection pts allocated to revision
  # who will receive either one-stage or two-stage revision based on the
  # surgeon's plan, i.e. these are not randomised.
  r_a_q <- data.table(
    silo = c("early", "late","late","late", "chronic", "chronic"),
    a = c("dair", "dair", "rev", "rev", "one", "two"),
    # intended/planned surgical approach
    qa = c("dair", "dair", "one", "two", "one", "two"),
    p = c(NA, NA, 0.5, 0.5, NA,  NA)
  )
  # eligibility for domain b
  r_eb <- data.table(
    qa = c("dair", "rev", "one", "two"),
    rand = c("N","Y","Y","Y")
  )
  r_b <- data.table(
    qa = c("dair", "one", "one", "two", "two"),
    b = c("w12", "w06p1", "w12p1", "d07p2", "w12p2")
  )
  r_b[, p := c(1, 0.5, 0.5, 0.5, 0.5)]
  r_ec <- data.table(
    silo = c("early", "late","chronic"),
    rand = c("Y","Y","Y"),
    p = rep(0.6, 3)
  )
  r_c <- data.table(
    c = c("norif", "rif"), p = rep(0.5, 2)
  )

  pkg.pop$spec <- list(
    r_silo = r_silo,
    r_joint = r_joint,
    r_ea = r_ea,
    r_a = r_a,
    r_a_q = r_a_q,
    r_eb = r_eb,
    r_b = r_b,
    r_ec = r_ec,
    r_c = r_c
  )

  pkg.pop$spec
}


# in case you want to override default spec
set_pop_spec <- function(spec) {
  pkg.pop$spec <- spec
}
