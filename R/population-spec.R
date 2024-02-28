library(data.table)

pkg.pop <- new.env()

get_pop_spec <- function(){

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

  pkg.pop$spec <- list(
    r_silo = r_silo,
    r_joint = r_joint,
    r_a = r_a,
    r_a_q = r_a_q,
    r_b = r_b,
    r_c = r_c
  )

  pkg.pop$spec
}




# in case you want to override default spec
set_pop_spec <- function(spec) {
  pkg.pop$spec <- spec
}
