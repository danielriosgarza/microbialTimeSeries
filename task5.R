##############Simulation conditions###################################

t.start = 0
t.end = 2000
t.step = 0.05
t.store = 500
norm = FALSE

simul.cond = list(t.start = t.start, t.end = t.end, t.step = t.step, t.store = t.store, norm = norm)

###############Common parameters#####################################
n.species = 25
migration.p = 0.1 
metacommunity.eveness = 1 #vary from 1 (uneven) to 1000 (highly even)
metacommunity.p = rdirichlet(1, rep(.5, n.species)*metacommunity.eveness)
growth.rates <- rbeta(n.species, shape1 = 0.5, shape2 = 0.5)
#growth.rates <- rep(1, n.species) 
x0 <- rep(0.001, n.species)
self.interactions <- runif(n.species, max(growth.rates)/3, max(growth.rates))

##############Stochastic ODE Models###################################
stochastic = FALSE
sigma.drift = 0.01
sigma.epoch = 0.1
sigma.external = 0.3
sigma.migration = 0.0
p.epoch = .1
n.external.perturbation = 1
t.external_events = getPerturbT(t.end, n.external.perturbation)
t.external_durations = rep(1, n.external.perturbation)

stoch.cond <- list(stochastic=stochastic, sigma.drift = sigma.drift,
                   sigma.epoch = sigma.epoch, sigma.external = sigma.external,
                   sigma.migration = sigma.migration, p.epoch = p.epoch, 
                   t.external_events = t.external_events,
                   t.external_durations = t.external_durations)


######################Params GLV Model ########################

####interaction matrix####
symmetric = FALSE
connectance = 0.05
diagonal = -1*self.interactions
mutualism.w = 1
commensalism.w = 1
parasitism.w = 1
amensalism.w = 1
competition.w = 1
interaction.w <- c(mutualism.w, commensalism.w, parasitism.w, amensalism.w,  competition.w)
scale = 0.5
sigma.A = 1
distribution = rnorm(n=n.species^2, sd = sigma.A)


rA_params <- list(n.species = n.species, connectance = connectance, 
                  diagonal = diagonal, interaction.w = interaction.w,
                  scale = scale, distribution = distribution, symmetric = symmetric)

A1 = eval(parse(text = "do.call(randomA, rA_params)"))

##########################

rA_params$connectance = 1
A2 = eval(parse(text = "do.call(randomA, rA_params)"))


A <- A1
#A[1,] = A2[1,]
#A[,1] = A2[,1]

A[1:5,] = A2[1:5,]
A[,1:5] = A2[,1:5]


#growth.rates[1] = 3*max(growth.rates)


#############################

glv_params <-  c(stoch.cond, simul.cond, list(n.species = n.species, A=A, 
                                              growth.rates = growth.rates,
                                              x0 = x0, migration.p=migration.p ))


glv <- parse(text = "do.call(simulateGLV, glv_params)")

######################Params Consumer Resource Model ########################

n.instances = 50
metaCom = matrix(0, n.instances, n.species)

for (i in 1:50){
  print(i)
  indices <- sample(1:n.species, 10)
  x0_ <- x0
  x0_[indices] = 0
  glv_params$x0 <- x0_
  simul <- eval(glv)
  metaCom[i,] <- simul$matrix[,colnames(simul$matrix)!="time"][t.store,]
  
  
}
gr <- as.numeric(metaCom[,1]==0) + 
   as.numeric(metaCom[,2]==0) + 
   as.numeric(metaCom[,3]==0) + 
   as.numeric(metaCom[,4]==0) +
   as.numeric(metaCom[,5]==0)

makeUMAP(metaCom, gradient = gr)
