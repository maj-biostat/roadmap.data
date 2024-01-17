\name{get_stan_data}
\alias{get_stan_data}
\title{Converts indexed data into list format}
\description{
  Creates data suitable for stan model currently used in sims.
}
\usage{
get_stan_data(d_i)
}
\details{
  Converts the indexed data to a list format suitable for the simulation stan model see \code{model-sim-02.stan} in sims.
}
\value{
  List returing aggregated data as \code{data.table} and corresponding stan data list comprising following entries.
  Indexed data is first aggregated and suitable for binomial likelihood.

  Early silo:

  \item{N_e  }{Number of rows of early silo data}
  \item{e_su }{Silo, joint infection index (1 to 6, see \code{get_indexes})}
  \item{e_y  }{Number of patients with treatment success}
  \item{e_n  }{Number of trials}
  \item{e_ec }{Membership indicator}
  \item{e_ecp}{Non-membership indicator}
  \item{e_c  }{Domain c treatment allocation}

  Late silo:

  \item{N_l  }{Number of rows of late silo data}
  \item{l_su }{Silo, joint infection index (1 to 6, see \code{get_indexes})}
  \item{l_y  }{Number of patients with treatment success}
  \item{l_n  }{Number of trials}
  \item{l_ec }{Membership indicator for domain c}
  \item{l_ecp}{Non-membership indicator for domain c}
  \item{l_c  }{Domain c treatment allocation}
  \item{l_ea  }{Membership indicator for domain a}
  \item{l_eap  }{Non-membership indicator for domain a}
  \item{l_a  }{Domain a treatment allocation}
  \item{l_eb1  }{Membership indicator for domain b (revision pt recving one-stage)}
  \item{l_eb2  }{Membership indicator for domain b (revision pt recving two-stage)}
  \item{l_ebp  }{Non-membership indicator for domain b}
  \item{l_b  }{Domain b treatment allocation}

  Chronic silo:

  \item{N_c  }{Number of rows of chronic silo data}
  \item{c_su }{Silo, joint infection index (1 to 6, see \code{get_indexes})}
  \item{c_y  }{Number of patients with treatment success}
  \item{c_n  }{Number of trials}
  \item{c_ec }{Membership indicator for domain c}
  \item{c_ecp}{Non-membership indicator for domain c}
  \item{c_c  }{Domain c treatment allocation}
  \item{c_ea  }{Membership indicator for domain a}
  \item{c_eap  }{Non-membership indicator for domain a}
  \item{c_a  }{Domain a treatment allocation}
  \item{c_eb1  }{Membership indicator for domain b (recving one-stage in domain a)}
  \item{c_eb2  }{Membership indicator for domain b (recving two-stage in domain a)}
  \item{c_ebp  }{Non-membership indicator for domain b}
  \item{c_b  }{Domain b treatment allocation}

}
