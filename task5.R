##############Simulation conditions###################################

t_start = 0
t_end = 2000
t_step = 0.05
t_store = 500
norm = FALSE

simul.cond = list(t_start = t_start, t_end = t_end, t_step = t_step, t_store = t_store, norm = norm)

###############Common parameters#####################################
n_species = 25
migration_p = 0.1 
metacommunity.eveness = 1 #vary from 1 (uneven) to 1000 (highly even)
metacommunity.p = rdirichlet(1, rep(.5, n_species)*metacommunity.eveness)
growth_rates <- rbeta(n_species, shape1 = 0.5, shape2 = 0.5)
#growth_rates <- rep(1, n_species) 
x0 <- rep(0.001, n_species)
self.interactions <- runif(n_species, max(growth_rates)/3, max(growth_rates))

##############Stochastic ODE Models###################################
stochastic = FALSE
sigma_drift = 0.01
sigma_epoch = 0.1
sigma_external = 0.3
sigma.migration = 0.0
p.epoch = .1
n_external.perturbation = 1
t_external_events = getPerturbT(t_end, n_external.perturbation)
t_external_durations = rep(1, n_external.perturbation)

stoch.cond <- list(stochastic=stochastic, sigma.drift = sigma.drift,
                   sigma_epoch = sigma_epoch, sigma_external = sigma_external,
                   sigma.migration = sigma.migration, p.epoch = p.epoch, 
                   t_external_events = t_external_events,
                   t_external_durations = t_external_durations)


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
distribution = rnorm(n=n_species^2, sd = sigma.A)


rA_params <- list(n_species = n_species, connectance = connectance, 
                  diagonal = diagonal, interaction.w = interaction.w,
                  scale = scale, distribution = distribution, symmetric = symmetric)

A1 = eval(parse(text = "do.call(randomA, rA_params)"))

##########################

rA_params$connectance = 1
A2 = eval(parse(text = "do.call(randomA, rA_params)"))


A <- A1
A[1,] = A2[1,]
A[,1] = A2[,1]

#A[1:5,] = A2[1:5,]
#A[,1:5] = A2[,1:5]


growth_rates[1:5] = 3*max(growth_rates)


#############################

glv_params <-  c(stoch.cond, simul.cond, list(n_species = n_species, A=A, 
                                              growth_rates = growth_rates,
                                              x0 = x0, migration_p=migration_p ))


glv <- parse(text = "do.call(simulateGLV, glv_params)")

######################Params Consumer Resource Model ########################

n.instances = 50
metaCom = matrix(0, n.instances, n_species)

for (i in 1:50){
  print(i)
  indices <- sample(1:n_species, 10)
  x0_ <- x0
  x0_[indices] = 0
  glv_params$x0 <- x0_
  simul <- eval(glv)
  metaCom[i,] <- simul$matrix[,colnames(simul$matrix)!="time"][t_store,]
  
  
}
gr <- as.numeric(metaCom[,1]==0) #+ 
#   as.numeric(metaCom[,2]==0) + 
#   as.numeric(metaCom[,3]==0) + 
#   as.numeric(metaCom[,4]==0) +
#   as.numeric(metaCom[,5]==0)

makeUMAP(metaCom, gradient = gr)
