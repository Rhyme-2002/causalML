<div align="center">

# causalML

### Causal Machine Learning Methods for Heterogeneous Treatment Effects

<p>
  <strong>Estimate • Simulate • Compare • Evaluate</strong>
</p>

<p>
  <a href="https://github.com/Rhyme-2002/causalML">
    <img src="https://img.shields.io/badge/GitHub-causalML-181717?style=for-the-badge&logo=github" alt="GitHub">
  </a>
  <img src="https://img.shields.io/badge/R-4.1%2B-276DC3?style=for-the-badge&logo=r" alt="R">
  <img src="https://img.shields.io/badge/License-MIT-green?style=for-the-badge" alt="License">
  <img src="https://img.shields.io/badge/Status-Under%20Development-orange?style=for-the-badge" alt="Status">
</p>

<p>
  An R package for estimating heterogeneous treatment effects using
  causal inference and machine learning methods.
</p>

</div>

---

## Overview

**causalML** is an R package designed for the estimation and evaluation
of heterogeneous treatment effects in observational data.

The package provides several causal inference and causal machine learning
methods for estimating:

- **Conditional Average Treatment Effect (CATE)**
- **Average Treatment Effect (ATE)**

It also provides a simulation framework for evaluating and comparing
different causal effect estimators under controlled data-generating
mechanisms.

---

## ✨ Features

| Feature | Description |
|:---|:---|
| 🎯 **CATE estimation** | Estimate heterogeneous treatment effects |
| 📊 **ATE estimation** | Estimate population-average treatment effects |
| 🌳 **Causal Forest** | Estimate heterogeneous effects using causal forests |
| 🌲 **Causal BART** | Bayesian Additive Regression Trees for causal effects |
| 🌿 **Causal Boosting** | Boosting based heterogeneous treatment effect estimation |
| 📐 **Causal MARS** | Multivariate Adaptive Regression Splines |
| 🌳 **Additive Causal Forest** | Additive treatment-effect estimation |
| 📈 **Parametric Standardization** | Logistic regression based standardization |
| 🧪 **Simulation** | Generate potential outcome probabilities |
| ⚖️ **Model Comparison** | Compare multiple causal estimators |
| 📏 **Performance Evaluation** | CATE RMSE and ATE ARB |
| 🔮 **Unified Prediction** | Common prediction interface for supported models |

---

# Installation

## Development Version

Install the development version directly from GitHub:

```r
install.packages("remotes")

remotes::install_github("Rhyme-2002/causalML")
