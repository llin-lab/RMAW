# RMAW

**Update: Nov 4 2024** 

A new and improved version of RMAW package will be further updated at this github page when necessary.

# Legacy version of RMAW

# RMAW - R package that computes the minimized aggregated Wasserstain distance and barycenter of Gaussian mixture models (GMMs). 

One key challenge encountered in single-cell data clustering is to combine clustering results of datasets acquired from multiple sources. We proposed to represent the clustering result of each data set by a Gaussian mixture model (GMM) and produced an integrated result based on the notion of Wasserstein barycenter. However, the precise barycenter of GMMs, a distribution on the same sample space, is computationally infeasible to solve. Importantly, the barycenter of GMMs may not be a GMM containing a reasonable number of components. We thus proposed to use the minimized aggregated Wasserstein (MAW) distance to approximate the Wasserstein metric and develop a new algorithm for computing the barycenter of GMMs under MAW. Recent theoretical advances further justify using the MAW distance as an approximation for the Wasserstein metric between GMMs. We also proved that the MAW barycenter of GMMs has the same expectation as the Wasserstein barycenter. Our proposed algorithm for clustering integration scales well with the data dimension and the number of mixture components, with complexity independent of data size. We demonstrated that the new method achieves better clustering results on several single-cell RNA-seq datasets than some other popular methods. We developed an R package to apply the algorithm above and compute MAW distance and barycenter.

For more informations please refer to the manuscript: [Lin Lin, Wei Shi, Jianbo Ye, Jia Li, Multisource Single-Cell Data Integration by MAW Barycenter for Gaussian Mixture Models, Biometrics, Volume 79, Issue 2, June 2023, Pages 866–877](https://academic.oup.com/biometrics/article/79/2/866/7513911)

# Install

Please refer to the github page to download and install Rmosek package before installing RMAW package: https://github.com/llin-lab/MAW

```R
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github('llin-lab/RMAW')
```

# Updates

**08.22.2024**: Modify README file to clarify installation steps.


**08.19.2024**: All the package functions are uploaded and all the dependency packages are included.


# Usage

```R
library(RMAW)

# Simplest use is running the wrapper function that creates a list object that contains the information of GMMs:

For more details on creating a list object see examples of each function.

# Contributors

RMAW was developed by Jingxuan Zhang and Lynn Lin. Please contact Jingxuan Zhang: jingxuan.zhang at duke edu for any questions or suggestions.

