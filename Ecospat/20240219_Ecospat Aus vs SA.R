################ 4. Ecological niche modelling - Ecospat ################
#### Packages and libraries ####


install.packages("ade4")
install.packages("devtools")
library(devtools) 
library(Rtools)

devtools::install_local("C:\\Users\\chris\\Desktop\\Meg full git repo\\20240220_ecospat\\ecospat\\ecospat\\ecospat", force = TRUE)
citation("ecospat")

library("ade4")
library("ecospat")
setwd('C:\\Users\\chris\\Desktop\\Meg full git repo\\20240220_ecospat')
#setwd("C:/Users/megan/OneDrive - Stellenbosch University/Chapter 2 - Tree ferns and climate change")

nat <- read.csv("20240220_occ_clim_Au.csv", sep = ",", dec = ".", header = TRUE)
inv <- read.csv("20240220_occ_clim_SA.csv", sep = ",", dec = ".", header = TRUE)




### PCA-ENVIRONMENT 
# the pca is calibrated on all the sites of the study area 
# Calibrating the PCA in the whole study area, including both native and invaded ranges (same as PCAenv in Broenniman et al. 2012)
pca.env <- dudi.pca(rbind(nat, inv)[,3:9], 
                    scannf = FALSE, 
                    nf = 2) 

# Plot Variables Contribution with ecospat.plot.contrib() - correlation circle
ecospat.plot.contrib(contrib = pca.env$co, 
                     eigen = pca.env$eig)

# predict the scores on the axes  
scores.globclim <- pca.env$li # PCA scores for the whole study area 

scores.sp.nat <- suprow(pca.env, nat[which(nat[,10]==1),3:9])$li # PCA scores for the species native distribution

scores.sp.inv <- suprow(pca.env, inv[which(inv[,10]==1),3:9])$li # PCA scores for the species invasive distribution 

scores.clim.nat <- suprow(pca.env, nat[,3:9])$li # PCA scores for the whole native study area 

scores.clim.inv <- suprow(pca.env, inv[,3:9])$li # PCA scores for the whole invaded study area


# Calculate the Occurrence Densities Grid with ecospat.grid.clim.dyn() 
# For a species in the native range (Australia) 
grid.clim.nat <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                       glob1 = scores.clim.nat, 
                                       sp = scores.sp.nat, 
                                       R = 100, #the resolution of the grid
                                       th.sp = 0) #the quantile used to delimit a threshold to exclude low species density values



# For a species in the invaded range (Reunion) 
grid.clim.inv <- ecospat.grid.clim.dyn(glob = scores.globclim, 
                                       glob1 = scores.clim.inv, 
                                       sp = scores.sp.inv, 
                                       R = 100, 
                                       th.sp = 0)


# Calculate Niche Overlap with ecospat.niche.overlap() 
niche_overlap <- ecospat.niche.overlap(grid.clim.nat, 
                      grid.clim.inv, 
                      cor = TRUE)$D 
print(niche_overlap)
#[1] 0.05324041

# Perform the Niche Equivalency Test with ecospat.niche.equivalency.test() 
eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, 
                                          grid.clim.inv, 
                                          rep = 100
                                          )
print(eq.test)

# Niche Similarity Test with ecospat.niche.similarity.test()
# shifting randomly the invasive niche in the invaded study area 

sim.test <- ecospat.niche.similarity.test(grid.clim.nat, 
                                          grid.clim.inv, 
                                          rep = 100, 
                                          intersection = 0.1) 
print(sim.test)

# Plot Equivalency test 
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")

# Plot Similarity test 
ecospat.plot.overlap.test(sim.test, "D", "Similarity")


# Delimiting niche categories and quantifying niche dynamics in analogue climates with ecospat.niche.dyn.index() 
niche_dynamics <- ecospat.niche.dyn.index(grid.clim.nat, 
                        grid.clim.inv, 
                        intersection = 0.1) 
print(niche_dynamics)
#result:
# $dynamic.index.w
# expansion stability unfilling 
# 0.0000000 1.0000000 0.7622818 


#Visualizing niche categories, niche dynamics and climate analogy between ranges with ecospat.plot.niche.dyn()
ecospat.plot.niche.dyn(grid.clim.nat, 
                       grid.clim.inv, 
                       quant = 0.1, 
                       interest = 2, 
                       title = "Niche Overlap", 
                       name.axis1 = "PC1", 
                       name.axis2 = "PC2")


# Plot Variables Contribution with ecospat.plot.contrib() 

ecospat.plot.contrib(contrib = pca.env$co, 
                     eigen = pca.env$eig)


# Plot the niche dynamics along one gradient (here bio16 Precipitation of Wettest Quarter) with ecospat.plot.niche.dyn()
# gridding the native niche 
grid.clim.t.nat <- ecospat.grid.clim.dyn(glob = as.data.frame(rbind(nat,inv)[,8]), 
                                         glob1 = as.data.frame(nat[,8]), 
                                         sp = as.data.frame(nat[which(nat[,10]==1),8]), 
                                         R = 100, 
                                         th.sp = 0)


# gridding the invaded niche 
grid.clim.t.inv <- ecospat.grid.clim.dyn(glob = as.data.frame(rbind(nat,inv)[,8]), 
                                         glob1 = as.data.frame(nat[,8]), 
                                         sp = as.data.frame(nat[which(nat[,10]==1),8]), 
                                         R = 100, 
                                         th.sp = 0)

t.dyn <- ecospat.niche.dyn.index(grid.clim.t.nat, 
                                 grid.clim.t.inv, 
                                 intersection=0.1) 

ecospat.plot.niche.dyn(grid.clim.t.nat, 
                       grid.clim.t.inv, 
                       quant = 0.1, 
                       interest = 2, 
                       title = "Niche Overlap", 
                       name.axis1 = "bio16 Precipitation of Wettest Quarter")



# Showing the shift of the niche centroid along the temperature gradient (compared to the shift of the available climate in the study area) with ecospat.shift.centroid() 
ecospat.shift.centroids(nat[which(nat[,10]==1),8], 
                        inv[which(inv[,10]==1),8],
                        nat[,8],
                        inv[,8])



