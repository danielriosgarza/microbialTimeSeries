---
title: "Consumer Resource Model - description"
author: "Daniel Garza"
date: "1/19/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## simulateConsumerResource

Function to simulate McArthur's consumer-resource model and its extensions for microbiomes. The model describes the dynamics of microbes and substrates in a community where the ecological interactions are explicitly encoded via the consumption and production of metabolites.

#### **Model Description**

Each microbial species from a community of $n$ different species has preferences for a specific subset of the $k$ different substrates that are available to the community. Their growth rates depend on the availability of these resources, through the following differential equation:

$\frac{dX_i}{dt} = \mu_{i}X_i (\sum_{j=1}^{k} e_{i,j} f_{i,j} S_j - \delta_i)$

where, where $\mu_{i}$ is the maximum growth rate of species $i$, while $X_i$ is its abundance at time $t$.

$e_{i,j}$ is the yield of substrate $j$ for species $i$, $f_{i,j}$ is the feeding form, which constraints the quantity of the substrate that can be consumed at each time step, while $S_j$ is the concentration of the substrate $j$ at time $t$, and $\delta_i$ is a dilution term.

The form of the feeding term is the Monod equation $f_{i,j}=\frac{S_j}{k_{i,j}+S_j}$, where $k_{i,j}$ is the Monod constant.

Resources change according to the following differential equation:

$\frac{dS_j}{dt}=\phi_j - \sum_{i=1}^{n} e_{i,j}f_{i,j}X_i + \sum_{i=1}^{n}p_{i,j} X_{i} (\sum_{j=1}^{k}e_{i,j}f_{i,j}-\delta_i)$

where $\phi_j$ is a dilution term and $p_{i,j}$ is the yield of production of substrate $j$ that results from the growth of species $i$.

In practice, a species either consumes or produces a metabolite (or is indifferent to its presence). Allowing us to summarize $e_{i,j}$ and $p_{i,j}$ in a single matrix ($E$, containing $n$ rows and $k$ columns). To distinguish them, the entries of $E$ are, respectively, positive and negative for the consumed and produced substrates.

There are many possibilities for structuring $E$ according to specific assumptions. We will later summarize the options that are built into miaSim. But, first let's jump to some quick examples of how to simulate a consumer resource model with miaSim.

### **Examples**

##### Example1: Default parameters

To illustrate the basic parameters, we begin simulating five species consuming or producing five different substrates. Some of the examples below are compatible with the Shiny app, click on the numbered example button on the top of the app screen to explore them.

If not provided by the user, all parameters have default values except for the the number of species and the number of resources, which need to be provided. Check the list of parameters below for a complete description of the parameters and their defaults

```{r}

n.species <- 5
n.resources <- 5

#simulate the model 

CRMsimul <- simulateConsumerResource(n.species = n.species, n.resources = n.resources)

#visualize the result 

makePlot(CRMsimul$matrix) #species plot
makePlotRes(CRMsimul$resources) #resources plot

```

##### Example2: exploring additional parameters

MiaSim provides a helper function to generate the matrix $E$.

Below is a simulation of the same model where the user has more control over the parameters.

```{r}

#generate the matrix E
E = randomE(n.species = n.species, n.resources = n.resources)
print(E)

#positive entries are consumed, negative produced, zero has no influence on the species

#define some simulation parameters
t.end = 2000 #when to stop
t.store = 500 #how many samples of the simulation to store (evenly spaced)
migration.p = 0 #whether to allow migration from a metacommunity
stochastic = 0 # whether to use noise in the simulation
dilution.rate = 0.001 #adding a dilution rate

#simulate the model
CRMsimul <- simulateConsumerResource(n.species = n.species, n.resources = n.resources, stochastic = stochastic, migration.p = migration.p, E=E, t.end = t.end, t.store = t.store, dilution.rate=dilution.rate)

#visualize the result
makePlot(CRMsimul$matrix) #species plot
makePlotRes(CRMsimul$resources) #resources plot

```

##### More about the $E$ matrix

The $E$ matrix contains the energy yield for the production of biomass and the secretion of metabolic by-products. In practice, there are alternative formulations to define the flux between substrate consumption and production.

When generating the random $E$ matrix using the function `randomE` there is a parameter called `maintenance` where the user may constraint the fraction of the flux that is not returned as a byproduct. This parameter can be interpreted as a fraction of the fluxes that are channeled into the organism's maintenance.

##### Example 3: maintenance

```{r}

E = randomE(n.species = 1, n.resources = 10, maintenance = .1)

print(sum(E*(E>0))) #consumed. For simplicity, values are normalized to add to 1

print(abs(sum(E*(E<0)))) #produced

```

##### Using stoichiometric relatioships for the secretion of metabolic by-products

MiaSim is flexible for other relationships between the consumption of growth compounds and the secretion of byproducts since the user may provide their own $E$ matrix to the `simulateConsumerResource` function.

One approach is to tie the consumption and production of resources to the stoichiometry of biological reactions (see [Marsland III, et al. 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006793)) .

##### Example 4: anaerobic food web

For instance, consider a community consisting of three microbes:

-   One that converts 1 mol of glucose to 3 mols of acetate with a yield of 4.3 mols of ATP, that we'll name "homoacetogen".

-   One that converts 1 mol of glucose into 2 mols of lactate, with a yield of two mols of ATP, that we'll name "homofermenter";

-   One that converts 4 mols of lactate and into 3 mols of butyrate, with a yield of 1 ATP, that we'll name "butyrateProducer".

```{r}
#The stoichiometric matrix 
D = matrix(c(1, -3, 0, 0, 1, 0, -2, 0, 0, 0, 4, -3), nrow = 3, byrow = TRUE)
yields = c(4.3/4, 2/4, 1/4)
E = D*yields

#growth rates
grs <- c(2, 4.5, 2.6)

#initial species composition
x0 <- c(1, 2, 1)

#initial media composition
resources <- c(10, 0, 0, 0)

#simulate the model
CRMsimul <- simulateConsumerResource(n.species = 3, n.resources = 4, stochastic = 0, migration.p = 0.0, E=E,dilution.rate=0.005,resources = resources, names.species = c('homoacetogenic', 'homofermentative', 'butyrateProducer'), names.resources = c('glucose', 'acetate', 'lactate', 'butyrate'), x0=x0, t.end = 10000, growth.rates = grs)

#visualize the result
makePlot(CRMsimul$matrix) #species plot
makePlot(CRMsimul$resources) #resources plot
```

##### Example 5: Preferred substrates

In general, substrates are not equally preferred by the microbes in a community, miaSim provides the possibility of generating the random $E$ matrix with a bias towards some preferred resource. This leads to the enrichment of the microbial preferences towards more "valuable" resources. This can be done by specifying the `trophic.preferences` parameter in the `randomE` function:

```{r}

E = randomE(n.species = 10, n.resources = 10, mean.consumption = 3, mean.production = 1, maintenance = 0.5, trophic.preferences = list(c(5, 3, 1, 1, 1, 1, 1, 1, 1, 1)))

#visualize the matrix
makeHeatmap(E)

```

##### Example 6: hierarchical trophic preferences

In a complex community, feeding might be structured in a hierarchical way through trophic groups. Each group would prefer certain substrates and secrete by-products that are then preferred by the next subgroup. Evidence of such multi-level organization has been shown for the human gut microbiome ([Wong et al. 2019](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007524)).

The parameter `trophic.levels` or the `randomE` function allows users to simulate communities with multiple levels of substrate preferences, that cross-feed to each other in hierarchical fashion.

```{r}
E = randomE(n.species = 20, n.resources = 20, mean.consumption = 3, mean.production = 2, maintenance = 0.0, trophic.levels = c(7, 13))

#visualize the matrix
makeHeatmap(E)

#visualize the sum of trophic preferences by level
Ep1 <- E[0:7,]
Ep2 <- E[7:20,]

levels<- t(cbind(colSums(Ep1*(Ep1<0)), colSums(Ep2*(Ep2>0))))

makeHeatmap(levels)

```

##### Example 6: Nested tropic levels

The default implementation of MiaSim supports only the definition of linear levels.

One can easily add nested relationships using the `trophic.preferences` parameter. For instance, consider the interactions depicted in the following cartoon where arrows represent the direction of the secretion/consumption flux of metabolic by-products between species.

![interaction cartoon](https://raw.githubusercontent.com/danielriosgarza/microbialTimeSeries/reviewed/files/images/interactionCartoon.png)

##### Adding Noise and simulating measurement error

##### More advanced example

##### Exploring the relationship between diversity and environmental noise

## .
