\name{get_pop_spec}
\alias{get_pop_spec}
\title{Create and retrieve data Population Specification}
\description{
  Creates a default list of population control variables.
}
\usage{
get_pop_spec()
}
\details{
  Population (design matrix) control variables are created .
  These are used as one of the arguments for \code{get_trial_data}.
}
\value{
  List comprising the following entries.

  \item{r_silo}{A 3 row data table capturing the expected proportion of sample by silo.}
  \item{r_joint}{A data table capturing the expected proportion of sample with each site of infection by silo.}
  \item{r_ea}{Indicator of entry into domain \code{a} by silo.}
  \item{r_a}{Allocation probabilities for each silo treatment options for domain \code{a}}
  \item{r_a_q}{Captures allocation probabilities for planned procedure for late stage patients allocated to revision in domain \code{a}}
  \item{r_eb}{Indicator of entry into domain \code{b} by silo.}
  \item{r_b}{Allocation probabilities for each domain \code{b} conditional on surgical procedure}
  \item{r_ec}{Indicator of entry into domain \code{b} by silo with expected proportion entering domain by silo.}
  \item{r_a}{Allocation probabilities for each silo treatment options for domain \code{c}}

}
