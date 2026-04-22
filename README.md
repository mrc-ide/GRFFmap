
# GRFFmap

This project uses a Gaussian Process (GP)-based spatial-temporal model to estimate the prevalence of antimalarial resistance markers across Africa over time.

### Overview
Resistance marker prevalence is modeled as a latent GP on a 2D spatial domain evolving over discrete yearly time steps. Observed data are binomial counts (mutant alleles out of total sequenced samples) linked to the latent field via a logistic transformation.

### Key Components
#### Spatial representation 
The spatial covariance follows a squared exponential kernel, approximated using Random Fourier Features (RFF) to reduce computational cost from O(N³) to a tractable linear operation. The latent logit-prevalence field is expressed as a linear combination of these features.
#### Temporal evolution 
Temporal dynamics are modeled as a Gaussian random walk over the RFF coefficients, forming a linear state-space model. This keeps computation linear in the number of time points O(T).
#### Inference 
Because the binomial likelihood is non-Gaussian, Pólya–Gamma data augmentation is used to convert observations into pseudo-Gaussian form, enabling exact Kalman filtering and RTS smoothing within an EM algorithm.
#### Output 
Posterior prevalence surfaces are reconstructed via Monte Carlo sampling, yielding estimates, 95% credible intervals, and exceedance probabilities at each grid cell and year.

### Installation

    git clone https://github.com/IDEELResearch/GRFFmap.git
    cd GRFFmap
    open GRFFmap.Rproj

`devtools::install_dev_deps()` will install all required packages, as
specified in the Imports in DESCRIPTION. (At a later date when analysis
is finalised, `renv` can be used to create a reproducible R environment
that anyone can use by calling `renv::restore` to set up package
dependencies.)

### Overview
The structure within analysis is as follows:

    R/                            # Packaged R functions 

    analysis/
        |
        ├── 01_xxxxx /            # analysis scripts used for generating figures
        |
        ├── data/               # data inputs for the model (shape files and prevalence data)

- Analysis scripts are to be run in the numbered order they are
  included. If there are shared numbers, then any order of these scripts
  works.

- Data that is **read only**, e.g. data shared from elsewhere and not
  generated using code in this repository, is stored in `data`

### Compendium DOI:

### Licenses

Code: [MIT](http://opensource.org/licenses/MIT) year: 2024, copyright
holder: OJ Watson

Data: [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse
