\name{get_sim_spec}
\alias{get_sim_spec}
\title{Create and retrieve Simulation Specification}
\description{
  Creates a default list of simulation control variables.
}
\usage{
get_sim_spec()
}
\details{
  Simulation control variables are created for each of the current model linear predictor components.
  These are used as one of the arguments for \code{get_trial_data}.
}
\value{
  List comprising the following entries.

  \item{a_s_u}{A named 3x2 matrix of baseline log-odds of response. Rows correspond to silos, columns to site of infection.}
  \item{gamma}{\code{gamm} is probably more complicated than it needs to be. A named 3x2x3 array of non-membership parameters. The last dimension of the array contains silos. Within each 3x2 matrix for each silo are domain (rows) and membership (Y/N). Typically Y entries will be set to zero and N (not randomised) will be set to some negative value. For the early silo only the c domain should be set. For the late and chronic silo, only the b and c domain should be set. In practice (but not in the simulations), non-membership would likely arise for late and chronic, but domain a has no randomisation options for domains a and b and therefore these should always be set to zero. Non-membership parameters are estimated as pooled across domains so the same value should be provided for each, e.g. domain c, not randomised could be set to -0.1 for all silos.}
  \item{b_a_late}{Parameters for late silo, domain a treatment effects with \code{dair} as reference group.}
  \item{b_a_chronic}{Parameters for chronic silo, domain a treatment effects with \code{one} as reference group.}
  \item{b_b1_late_one}{Parameters for late silo, domain b treatment effects for patients randomised to revision and with planned one-stage procedure in domain a with \code{w06p1} as reference group.}
  \item{b_b1_late_two}{Parameters for late silo, domain b treatment effects for patients randomised to revision and with planned two-stage procedure in domain a with \code{d07p2} as reference group.}
  \item{b_b1_chronic_one}{Parameters for chronic silo, domain b treatment effects for patients randomised to one-stage procedure in domain a with \code{w06p1} as reference group.}
  \item{b_b2_chronic_two}{Parameters for chronic silo, domain b treatment effects for patients randomised to two-stage procedure in domain a with \code{d07p2} as reference group.}
  \item{b_c}{Parameters (shared across all silos) for domain c treatment effects in domain c with \code{norif} as reference group.}

}
