# Package `spyce`
## How to install
*  Easiest: install from source on GitHub
```{r}
install.packages("remotes")
remotes::install_github("khnhan/spyce")
```

* If you see "HTTP error 401. Bad credentials":
```{r}
install.packages("gitcreds")
gitcreds::gitcreds_delete()  # clears a bad cached token
remotes::install_github("khnhan/spyce")
```
## Description
This package analyzes survival data with outcome-dependent right-censored covariate. The semiparametric estimator SPYCE is implemented via `get_SPYCE_...` functions. This package provides other estimators for analyzing the same data: complete case estimator (CC), imputation estimator, inverse probability weighting estimator (IPW), and the maximum likelihood estimator (MLE). The estimated asymptotic variance for the estimators can be computed with corresponding `variance_...` functions. Conditional Kaplan-Meier estimator as an auxiliary tool for nonparametric estimators can also be implemented with this package.

## Usage
The full data is $(Y,X,C,Z)$, and the observed data is $(Y,W,\Delta,Z)$, where $W=\min(X,C)$ and $\Delta = 1(X\le C)$.
It is assumed that $X\bot C\mid Y,Z$, and $\eta_1 = f_{X|Z}$ and $\eta_2 = f_{C|Y,Z}$ are considered nuisance distributions.

The functions `get_...` find an estimator of $\beta$ and $\sigma$, and 
the functions `variance_...` compute the estimated asymptotic variance of the estimators.:
- With `z_data` present, then $Y\mid X,Z \sim N\big(\beta_0+\beta_1 X+\beta_2 Z+\beta_3 ZX,\ \sigma^2\big)$.
- Without `z_data`, then $Y\mid X \sim N\!\big(\beta_0+\beta_1 X,\ \sigma^2\big)$.

Each of the estimators (CC, imputation, IPW, MLE, SPYCE) has a parametric nuisance distribution version (truncated normal) and a nonparametric nuisance distribution version, while each supports four modes:
- (A) estimate $(\beta,\sigma)$ with `z_data`
- (B) estimate $(\beta,\sigma)$ without `z_data`
- (C) estimate $\beta$ with `z_data` and fixed $\sigma$
- (D) estimate $\beta$ without `z_data` and fixed $\sigma$

For the more detailed example, see **Examples** in `help(get_SPYCE_param12)` and `help(variance_SPYCE_param12)`.