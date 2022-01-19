---
title: "Consumer Resource Model"
author: "Daniel Garza"
date: "1/19/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## simulateConsumerResource

Function to simulate McArthur's consumer-resource model and extensions for microbiomes. The model describes the explicit interaction between microbes in a community via the consumption and production of metabolites.

#### **Model Description**

Each microbial species from a pool of $n$ different species has preferences for specific subset of the $k$ different substrates that are available to the community. Their growth rates depend on the availability of these resources, through the following differential equation:

$\frac{dX_i}{dt} = \mu_{i}X_i (\sum_{j=1}^{k} e_{i,j} f_{i,j} S_j - \delta_i)$

where, where $\mu_{i}$ is the maximum growth rate of species $i$, while $X_i$ is its abundance at time $t$.

$e_{i,j}$ is the yield of substrate $j$ for species $i$, $f_{i,j}$ is the feeding form, which constraints the quantity of the substrate that can be consumed at each time step, while $S_j$ is the concentration of the substrate $j$ at time $t$, and $\delta_i$ is a dilution term.

In the current implementation the yields are be provided as a matrix $E$ (containing $n$ rows and $k$ columns). The form of the feeding term is the Monod equation $f_{i,j}=\frac{S_j}{k_{i,j}+S_j}$, where $k_{i,j}$ is the Monod constant.

Resources change according to the following differential equation:

$\frac{dS_j}{dt}=\phi_j - \sum_{i=1}^{n} e_{i,j}f_{i,j}X_i + \sum_{i=1}^{n}p_{i,j} X_{i} (\sum_{j=1}^{k}e_{i,j}f_{i,j}-\delta_i)$

where $\phi_j$ is a dilution term and $p_{i,j}$ is the yield of production of substrate $j$ that results from the growth species $i$.

In practice, a species either consumes or produces a metabolite (or is indifferent to its presence). Allowing us to summarize $e_{i,j}$ and $p_{i,j}$ in a single matrix ($E$). To distinguish them, the entries of $E$ are, respectively, positive and negative for the consumed and produced substrates.

There are many possibilities of for structuring $E$ according to specific assumptions. We will later summarize the options that are built into miaSim. We will first quickly provides

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
