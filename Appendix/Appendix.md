# Appendix

## Conditional Posterior Distribution

The conditional posterior distributions for the common parameters are given as follows:

$$
\begin{aligned}
\pi(\theta_{k})\propto & \left[\prod_{K=1}^{N}\prod_{i=1}^{P}\mathbb{P}(Y_{k,i}=y_{k,i}|\boldsymbol{\Theta})\right]\times
  \left[N_{\theta_{k}}(0,\sigma^{2})\right]\\
\pi(\beta_{i})\propto & \left[\prod_{k=1}^{N}\prod_{i=1}^{P}\mathbb{P}(Y_{k,i}=y_{k,i}|\boldsymbol{\Theta})\right]\times  \left[N_{\beta_{i}}(0,\tau_{\beta}^{2})\right]\\  
\pi(\gamma) \propto & \left[\prod_{k=1}^{N}\prod_{i=1}^{P}\mathbb{P}(Y_{k,i}=y_{k,i}|\boldsymbol{\Theta})\right]\times \left[\text{Log-normal}_{\gamma}(\mu_{\gamma},\tau_{\gamma}^2)\right]\\
\pi(\boldsymbol{z_{k}})\propto & \left[\prod_{k=1}^{N}\prod_{i=1}^{P}\mathbb{P}(Y_{k,i}=y_{k,i}|\boldsymbol{\Theta})\right]\times
  \left[\text{MVN}_{D,\boldsymbol{z_{k}}}(\boldsymbol{0,I_{D}})\right]\\
\pi(\boldsymbol{w_{i}})\propto & \left[\prod_{k=1}^{N}\prod_{i=1}^{P}\mathbb{P}(Y_{k,i}=y_{k,i}|\boldsymbol{\Theta})\right]\times
  \left[\text{MVN}_{D,\boldsymbol{w_{i}}}(\boldsymbol{0,I_{D}})\right]\\
\pi(\sigma^{2})
\propto &\,\text{Inv-Gamma}\left(\left(\frac{N}{2}+a_{\sigma}\right),\frac{1}{2}\sum_{k=1}^{N}\theta_{k}^{2}+b_{\sigma}\right).
\end{aligned}
$$

The conditional posterior distributions of $\alpha$ and $\sigma_{\epsilon}^2$ for the 2PL LSIRM and the LSIRM-continuous, respectively, are given as follows:

$$
\begin{aligned}
    \pi(\alpha_{i})\propto & \left[\prod_{k=1}^{N}\prod_{i=1}^{P}\mathbb{P}(Y_{k,i}=y_{k,i}|\theta_{k},\alpha_{i}, \beta_{i},\gamma,\boldsymbol{z_{k},w_{i}})\right]\times  \left[\text{Log-normal}_{\alpha_{i}}(0,\tau_{\alpha}^{2})\right]\\
    \pi(\sigma_{\epsilon}^{2})
    \propto & \,\text{Inv-Gamma}\left(\left(\frac{NP}{2}+a_{\sigma_\epsilon}\right),\frac{1}{2}\sum_{k=1}^{N}\sum_{i=1}^{P}\left(y_{k,i}-\left(\theta_k +\beta_i - \gamma || {\boldsymbol{z_k}} - {\boldsymbol{w_i}} ||\right)\right)^{2}+b_{\sigma_\epsilon}\right),
\end{aligned}
$$

where $\boldsymbol{\Theta} = \{\theta_{k},\beta_{i},\gamma,\boldsymbol{z_{k},w_{i}}\}$ for 1PL LSIRM, $\boldsymbol{\Theta} = \{\theta_{k}, \alpha_{i}, \beta_{i}, \gamma, \boldsymbol{z_{k},w_{i}}\}$ for 2PL LSIRM, and $\boldsymbol{\Theta} = \{\theta_{k}, \beta_{i}, \gamma, \boldsymbol{z_{k},w_{i}}, \sigma_{\epsilon}^2\}$ for 1PL LSIRM-continuous.

## Default Value for the Prior Specification and the Jumping Rules

| Parameter | Models | Arguments | Default value |
|-----------|--------|-----------|---------------|
| $\theta$  | `lsirm1pl`, `lsirm2pl` | `pr_mean_theta` | 0 |
| $\beta$   | `lsirm1pl`, `lsirm2pl` | `pr_mean_beta`, `pr_sd_beta` | 0, 1 |
| $\alpha$  | `lsirm2pl` | `pr_mean_alpha`, `pr_sd_alpha` | 0.5, 1 |
| $\log\gamma$ | `lsirm1pl`, `lsirm2pl` | `pr_mean_gamma`, `pr_sd_gamma` | 0.5, 1 |
| $\sigma^{2}$ | `lsirm1pl`, `lsirm2pl` | `pr_a_theta`, `pr_b_theta` | 0.001, 0.001 |
| $\sigma_{\epsilon}^{2}$ | `lsirm1pl`, `lsirm2pl` | `pr_a_eps`, `pr_b_eps` | 0.001, 0.001 |

*Table 1: Model parameters, related models, corresponding arguments, and default values for priors. Note that for $\log \gamma$, the mode of its prior distribution is 0.61, with mean 2.72, and standard deviation 3.56.*

| Parameter | Models | Arguments | Default value |
|-----------|--------|-----------|---------------|
| $\theta$  | `lsirm1pl`, `lsirm2pl` | `jump_theta` | 1 |
| $\beta$   | `lsirm1pl`, `lsirm2pl` | `jump_beta` | 0.4 |
| $\alpha$  | `lsirm2pl` | `jump_alpha` | 1 |
| $\log\gamma$ | `lsirm1pl`, `lsirm2pl` | `jump_gamma` | 0.025 |
| $\boldsymbol{z}$ | `lsirm1pl`, `lsirm2pl` | `jump_z` | 0.5 |
| $\boldsymbol{w}$ | `lsirm1pl`, `lsirm2pl` | `jump_w` | 0.5 |

*Table 2: Model parameters, related models, corresponding arguments, and default values for jumping rules.*

## List of Package Argument

| Arguments | Description |
|-----------|-------------|
| `data` | Item response matrix to be analyzed. |
| `ndim` | Dimension of latent space. |
| `niter` | Number of iterations to run MCMC sampling. |
| `nburn` | Number of initial, pre-thinning, MCMC iterations to discard. |
| `nthin` | Number of thinning, MCMC iterations to discard. |
| `nprint` | MCMC samples are displayed. |
| `pr_mean_theta` | Mean of the normal prior distribution for $\boldsymbol{\theta}$. |
| `pr_mean_beta` | Mean of the normal prior distribution for $\boldsymbol{\beta}$. |
| `pr_sd_beta` | Standard deviation of the normal prior distribution for $\boldsymbol{\beta}$. |
| `pr_mean_alpha` | Mean of the normal prior distribution for $\alpha$. |
| `pr_sd_alpha` | Standard deviation of the normal prior distribution for $\alpha$. |
| `pr_mean_gamma` | Mean of the log-normal prior distribution for $\gamma$. |
| `pr_sd_gamma` | Standard deviation of the log-normal prior distribution for $\gamma$. |
| `pr_a_theta` | Shape parameter of the inverse gamma prior for the variance of $\boldsymbol{\theta}$, $\sigma$. |
| `pr_b_theta` | Scale parameter of the inverse gamma prior for the variance of $\boldsymbol{\theta}$, $\sigma$. |
| `pr_a_eps` | Shape parameter of the inverse gamma prior for the variance of data likelihood, $\sigma_{\epsilon}^2$. |
| `pr_b_eps` | Scale parameter of the inverse gamma prior for the variance of data likelihood, $\sigma_{\epsilon}^2$. |
| `jump_theta` | Jumping rule of the proposal density for $\boldsymbol{\theta}$. |
| `jump_beta` | Jumping rule of the proposal density for $\boldsymbol{\beta}$. |
| `jump_alpha` | Jumping rule of the proposal density for $\alpha$. |
| `jump_gamma` | Jumping rule of the proposal density for $\gamma`. |
| `jump_z` | Jumping rule of the proposal density for $z$. |
| `jump_w` | Jumping rule of the proposal density for $w`. |
| `missing` | Replaced value for missing data. |

*Table 3: Arguments for prior distributions, jumping rules, and others.*

| Arguments | Description |
|-----------|-------------|
| `data` | Item response matrix to be analyzed. |
| `bic` | Bayesian information criterion value. |
| `mcmc_inf` | Number of MCMC iteration, burn-in periods, and thinning intervals. |
| `map_inf` | Log maximum a posteriori probability estimates and acceptance rates for all parameters. |
| `post_means` | Posterior means for all parameters. |
| `post_sds` | Posterior standard deviation for all parameters. |
| `theta` | Estimates of the person parameters. |
| `beta` | Estimates of the item intercept parameters. |
| `alpha` | Estimates of the item slope parameters (only for 2PL LSIRM). |
| `gamma` | Estimates of the variance of latent positions. |
| `z` | Estimates of the latent positions of persons. |
| `w` | Estimates of the latent positions of items. |

*Table 4: Posterior inference and model fit information.*
