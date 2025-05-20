%% File Name: mnlfa.Rd
%% File Version: 0.380

\name{mnlfa}
\alias{mnlfa}
\alias{summary.mnlfa}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Moderated Nonlinear Factor Analysis
}

\description{
General function for conducting moderated nonlinear factor analysis
(Curran et al., 2014). Item slopes and item intercepts can be modeled as functions
of person covariates.

Parameter regularization is allowed. For categorical covariates,
group lasso can be used for regularization.
}


\usage{
mnlfa(dat, items, weights=NULL, item_type="2PL", formula_int=~1, formula_slo=~1,
   formula_res=~0, formula_mean=~0, formula_sd=~0, theta=NULL, parm_list_init=NULL,
   parm_trait_init=NULL, prior_init=NULL, regular_lam=c(0,0,0), regular_alpha=c(0,0,0),
   regular_type=c("none", "none", "none"), maxit=1000, msteps=4, conv=1e-05,
   conv_mstep=1e-04, h=1e-04, parms_regular_types=NULL, parms_regular_lam=NULL,
   parms_regular_alpha=NULL, parms_iterations=NULL, center_parms=NULL, center_max_iter=6,
   L_max=.07, verbose=TRUE, numdiff=FALSE)

\method{summary}{mnlfa}(object, file=NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame with item responses
}
  \item{items}{
Vector containing item names
}
\item{weights}{Optional vector of sampling weights for persons}
  \item{item_type}{
String or vector of item types. The item types \code{"1PL"} or
\code{"2PL"} for dichotomous items, \code{"GPCM"} for polytomous and \code{"NO"} for
continuous items can be chosen.
}
  \item{formula_int}{
String or list with formula for item intercepts
}
  \item{formula_slo}{
String or list with formula for item slopes
}
  \item{formula_res}{
String or list with formula for logarithms of residual standard deviations
}
  \item{formula_mean}{
Formula for mean of the trait distribution
}
  \item{formula_sd}{
Formula for standard deviation of the trait distribution
}
  \item{theta}{
Grid of \eqn{\theta} values used for approximation of normally
distributed trait
}
  \item{parm_list_init}{
Optional list of initial item parameters
}
\item{parm_trait_init}{Optional list of initial parameters for trait distribution}
\item{prior_init}{Optional matrix of prior distribution for persons}
  \item{regular_lam}{
Vector of length 2 or 3 (for \code{item_type="NO"}) containing regularization
parameters \eqn{\lambda} for item intercepts and item slopes
}
\item{regular_alpha}{Vector of length 2 or 3 containing \eqn{\alpha}
regularization parameter}
  \item{regular_type}{
Type of regularization method. Can be \code{"none"},
\code{"lasso"}, \code{"scad"}, \code{"mcp"}, \code{"scadL2"}, \code{"ridge"} or
\code{"elnet"}.
}
  \item{maxit}{
Maximum number of iterations
}
  \item{msteps}{
Maximum number of M-steps
}
  \item{conv}{
Convergence criterion with respect to parameters
}
  \item{conv_mstep}{
Convergence criterion in M-step
}
  \item{h}{
Numerical differentiation parameter
}
  \item{parms_regular_types}{
Optional list containing parameter specific regularization types
}
  \item{parms_regular_lam}{
Optional list containing parameter specific regularization parameters
}
  \item{parms_regular_alpha}{
Optional list containing parameter specific regularization parameters
}
  \item{parms_iterations}{
Optional list containing sequence of parameter indices used for updating
}
\item{center_parms}{Optional list indicating which parameters
should be centered during initial iterations.}
\item{center_max_iter}{Maximum number of iterations in which parameters
  should be centered.}
\item{L_max}{Majorization parameter used in regularization}
\item{verbose}{
Logical indicating whether output should be printed
}
\item{numdiff}{Logical indicating whether numerical differentiation
should be used}
\item{object}{Object of class \code{mnlfa}}
\item{file}{Optional file name}
\item{\dots}{Further arguments to be passed}
}

\details{
The moderated factor analysis model for dichotomous responses defined as
\deqn{P(X_{pi}=1 | \theta_p )=invlogit( a_{pi} \theta_p - b_{pi} ) }
The trait distribution \eqn{\theta_p \sim N( \mu_p, \sigma_p^2)}
allows a latent regression of person covariates on the mean
with \eqn{\mu_p=\bold{X}_p \bold{\gamma}} (to be specified in \code{formula_mean})
and the logarithm of the standard deviation \eqn{\log \sigma_p=\bold{Z}_p \bold{\delta} }
(to be specified in \code{formula_sd}).
Item intercepts and item slopes can be moderated by person covariates, i.e.
\eqn{a_{pi}=\bold{W}_{pi} \bold{\alpha}_i } and
\eqn{b_{pi}=\bold{V}_{pi} \bold{\beta}_i }. Regularization on (some of) the
\eqn{\bold{\alpha}_i} or \eqn{\bold{\beta}_i} parameters is allowed.

For polytomous item responses, the generalized partial credit model is parametrized
as \deqn{P(X_{pi}=k | \theta_p \propto \exp ( k a_{pi} \theta_p - k b_{pi} - b_{k} }
with \eqn{b_0=0}.

For normally distributed responses, the conditional distribution of item responses
is defined as \deqn{ X_{pi} | \theta_p \sim \mathrm{N} ( b_{pi} + a_{pi} \theta_p,
\psi_{pi}^2 ) } Note that \eqn{\log \psi_{pi}} is modeled in this function.

The model is estimated using an EM algorithm with the coordinate descent
method during the M-step (Sun et al., 2016).
}


\value{
List with model results including

\item{item}{Summary table for item parameters}
\item{trait}{Summary table for trait parameters}
}


\references{
Curran, P. J., McGinley, J. S., Bauer, D. J., Hussong, A. M., Burns, A.,
Chassin, L., Sher, K., & Zucker, R. (2014).
A moderated nonlinear factor model for the development of commensurate measures
in integrative data analysis. \emph{Multivariate Behavioral Research, 49}(3),
214-231.
http://dx.doi.org/10.1080/00273171.2014.889594

Sun, J., Chen, Y., Liu, J., Ying, Z., & Xin, T. (2016). Latent variable selection for
multidimensional item response theory models via L1 regularization.
\emph{Psychometrika, 81}(4), 921-939.
 https://doi.org/10.1007/s11336-016-9529-6
}


%\author{
%%  ~~who you are~~
%}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also the \pkg{aMNLFA} package for automatized moderated nonlinear factor analysis
which provides convenient wrapper functions for automized analysis
in the \emph{Mplus} software.

See the \pkg{GPCMlasso} package for the regularized generalized partial credit model.
}

\examples{
#############################################################################
# EXAMPLE 1: Dichotomous data, 1PL model
#############################################################################

data(data.mnlfa01, package="mnlfa")

dat <- data.mnlfa01
# extract items from dataset
items <- grep("I[0-9]", colnames(dat), value=TRUE)
I <- length(items)

# maximum number of iterations (use only few iterations for the only purpose of
# providing CRAN checks)
maxit <- 10

#***** Model 1: 1PL model without moderating parameters and without covariates for traits

# no covariates for trait
formula_mean <- ~0
formula_sd <- ~1
# no item covariates
formula_int <- ~1
formula_slo <- ~1

mod1 <- mnlfa::mnlfa( dat=dat, items, item_type="1PL", formula_int=formula_int,
             formula_slo=formula_slo, formula_mean=formula_mean, formula_sd=formula_sd,
             maxit=maxit )
summary(mod1)

\donttest{
#***** Model 2: 1PL model without moderating parameters and with covariates for traits

# covariates for trait
formula_mean <- ~female + age
formula_sd <- ~1

mod2 <- mnlfa::mnlfa( dat=dat, items, item_type="1PL", formula_int=formula_int,
             formula_slo=formula_slo, formula_mean=formula_mean, formula_sd=formula_sd)
summary(mod2)

#***** Model 3: 1PL model with moderating parameters and with covariates for traits
#***   Regularization method 'mcp'

# covariates for trait
formula_mean <- ~female + age
formula_sd <- ~1
# moderation effects for items
formula_int <- ~1+female+age
formula_slo <- ~1

# center parameters for female and age in initial iterations for improving convergence
center_parms <- list( rep(2,I), rep(3,I) )

# regularization parameters for item intercept and item slope, respectively
regular_lam <- c(.06, .25)
regular_type <- c("mcp","none")

mod3 <- mnlfa::mnlfa( dat=dat, items, item_type="1PL", formula_int=formula_int,
            formula_slo=formula_slo, formula_mean=formula_mean, formula_sd=formula_sd,
            center_parms=center_parms, regular_lam=regular_lam, regular_type=regular_type )
summary(mod3)

# update Model 3 with initial parameters
parm_trait_init <- mod3$parm_trait
parm_list_init <- mod3$parm_list
prior_init <- mod3$prior
mod3b <- mnlfa::mnlfa( dat=dat, items, item_type="1PL", formula_int=formula_int,
            formula_slo=formula_slo, formula_mean=formula_mean, formula_sd=formula_sd,
            center_parms=center_parms, regular_lam=regular_lam,
            regular_type=regular_type, parm_trait_init=parm_trait_init,
            parm_list_init=parm_list_init, prior_init=prior_init,  )
summary(mod3b)

#***** Model 4: 1PL model with selected moderated item parameters

#* trait distribution
formula_mean <- ~0+female+age
formula_sd <- ~1

#* formulas for item intercepts
formula_int <- ~1
formula_int <- mnlfa::mnlfa_expand_to_list(x=formula_int, names_list=items)
mod_items <- c(4,5,6,7)
for (ii in mod_items){
    formula_int[[ii]] <- ~1+female+age
}
formula_slo <- ~1

mod4 <- mnlfa::mnlfa( dat=dat, items, item_type="1PL", formula_int=formula_int,
              formula_slo=formula_slo, formula_mean=formula_mean, formula_sd=formula_sd)
mod4$item
mod4$trait
summary(mod4)

#############################################################################
# EXAMPLE 2: Continuous items
#############################################################################

#-- simulate data
set.seed(78)
N <- 1000
I <- 10
z <- stats::runif(N, -1, 1)
theta <- stats::rnorm(N, mean=0.3*z - .1*z^2, sd=exp(-.5+.3*z) )
dat <- matrix(NA, nrow=N, ncol=I)
items <- colnames(dat) <- paste0("I",1:I)
dat <- as.data.frame(dat)
dat$z <- z

for (ii in 1L:I){
    f <- 0.4*(ii==1)
    dat[,ii] <- .6+-f*z^2 + (1+f*z)*theta + stats::rnorm(N, sd=exp(-.3)  )
}

#--- estimate model
formula_mean <- ~ 0 + z
formula_sd <- ~ 0 + z
formula_int <- ~1+z
formula_slo <- ~1+z
formula_res <- ~1+z

regular_type <- c("mcp","mcp","mcp")
regular_lam <- rep(1e-3,3)
maxit <- 10
item_type <- "NO"

mod1 <- mnlfa::mnlfa( dat=dat, items, item_type=item_type, formula_int=formula_int,
             formula_slo=formula_slo, formula_res=formula_res, formula_mean=formula_mean,
             formula_sd=formula_sd, maxit=maxit, regular_type=regular_type, h=1e-4,
             regular_lam=regular_lam)
summary(mod1)

#############################################################################
# EXAMPLE 3: Polytomous items
#############################################################################

#--- simulate data
set.seed(78)
N <- 2000
I <- 8
z <- stats::runif(N, -1, 1)
theta <- stats::rnorm(N, mean=0.3*z, sd=exp(-0.3*z) )
dat <- matrix(NA, nrow=N, ncol=I)
items <- colnames(dat) <- paste0("I",1:I)
dat <- as.data.frame(dat)
dat$z <- z
for (ii in 1L:I){
    f <- 0.4*(ii \%in\% c(1,3) )
    y <- -f*z + (1+f*z)*theta + stats::rnorm(N, sd=exp(-.3)  )
    K <- 2 + (ii \%\% 2==0 )
    br <- c(-Inf,seq(-1,1, len=K), Inf) + ii / I
    y <- as.numeric( cut( y, breaks=br) )-1
    dat[,ii] <- y
}

#--- estimate model
formula_mean <- ~1 + z
formula_sd <- ~1 + z
formula_int <- ~1+z
formula_slo <- ~1+z

regular_type <- rep("scadL2",2)
regular_lam <- rep(0.02,2)
regular_alpha <- rep(0.5,2)
maxit <- 10
item_type <- "GPCM"

mod1 <- mnlfa::mnlfa( dat=dat, items, item_type=item_type, formula_int=formula_int,
             formula_slo=formula_slo, formula_mean=formula_mean,
             formula_sd=formula_sd, maxit=maxit, regular_type=regular_type, h=1e-4,
             regular_lam=regular_lam, regular_alpha=regular_alpha)
summary(mod1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
