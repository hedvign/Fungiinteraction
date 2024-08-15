# Fungiinteraction
Data and models for identifying fungi interactions

Combining observational and experimental data to estimate environmental and species drivers of fungal metacommunity dynamics (in revision)

Hedvig Nenzen


This data was previously published, see these sources for more details: 

Jönsson, M., Edman, M. & Jonsson, B.G. (2008) Colonization and extinction patterns of wood-decaying fungi in a boreal old-growth Picea abies forest. Journal of Ecology, 96, 1065–1075.

Ovaskainen, O., Hottola, J. & Shtonen, J. (2010) Modeling species co-occurrence by multivariate logistic regression generates new hypotheses on fungal interactions. Ecology, 91, 2514–2521.

Ottosson, E., Nordén, J., Dahlberg, A., Edman, M., Jönsson, M., Larsson, K.H., Olsson, J., Penttilä, R., Stenlid, J. & Ovaskainen, O. (2014) Species associations during the succession of wood-inhabiting fungal communities. Fungal Ecology, 11, 17–28.


Data Collection Overview:

We pooled data from two studies that surveyed the presence or absence of wood-decay fruit bodies on a total of 1379 dead Norway spruce (Picea abies) logs. The logs were surveyed twice, 6 years apart. The study sites were composed of old-growth unmanaged forests where we studied naturally fallen logs. In Brattiken nature reserve (Sweden, 65°25"N, 16°06´E), 843 logs were inventoried in autumn of 1997 and 2003 (Edman & Jonsson, 2001; Berglund et al., 2005; Jönsson et al., 2008). In Rörstrand nature reserve (Finland, 60°27"N, 25°11´E), 536 logs were inventoried in autumn of 2002 and 2008 (Ovaskainen et al., 2010; Ottosson et al., 2014). Among a total of 210 fungal species, we chose 12 species that were present at both study sites in sufficient numbers and were expected to interact based on the literature. We used additional repeat visit data to account for imperfect detection (as Moor et al., 2021). Dead fruit bodies were excluded from the analysis. We included two log-level variables that are well-known to explain species occurrences: log decay stage and log diameter (cm). All environmental variables were transformed to mean zero and unit standard deviation (subtracting the mean and dividing each value by the standard deviation). We also included the decay stage variable squared, as some species are most frequent (in terms of occupancy, i.e.  number of presences divided by number of logs) in logs of intermediate decay.


Species: 
Fomitopsis pinicola, 
Heterobasidion parviporum
Trichaptum abietinum
Phlebia centrifuga
Phellinidium ferrugineofuscum
Fomitopsis rosea
Postia caesia
Neoantrodia serialis
Fuscoporia viticola
Gloeophyllum sepiarium
Skeletocutis brevispora
Phellopilus nigrolimitatus

References:

Berglund, H., Edman, M. & Ericson, L. (2005) Temporal variation of wood-fungi diversity in boreal old-growth forests: Implications for monitoring. Ecological Applications, 15, 970–982.

Edman, M. & Jonsson, B.G. (2001)  Spatial pattern of downed logs and wood‐decaying fungi in an old‐growth Picea abies forest . Journal of Vegetation Science, 12, 609–620.

Jönsson, M., Edman, M. & Jonsson, B.G. (2008) Colonization and extinction patterns of wood-decaying fungi in a boreal old-growth Picea abies forest. Journal of Ecology, 96, 1065–1075.

Moor, H., Nordén, J., Penttilä, R., Siitonen, J. & Snäll, T. (2021) Long-term effects of colonization–extinction dynamics of generalist versus specialist wood-decaying fungi. Journal of Ecology, 109, 491–503.

Ovaskainen, O., Hottola, J. & Shtonen, J. (2010) Modeling species co-occurrence by multivariate logistic regression generates new hypotheses on fungal interactions. Ecology, 91, 2514–2521.

Ottosson, E., Nordén, J., Dahlberg, A., Edman, M., Jönsson, M., Larsson, K.H., Olsson, J., Penttilä, R., Stenlid, J. & Ovaskainen, O. (2014) Species associations during the succession of wood-inhabiting fungal communities. Fungal Ecology, 11, 17–28.


Model:
We developed metacommunity models where the response variable is log colonization or extinction of each individual species. To statistically explain each species’ dynamics, explanatory variables were environmental variables and the presence of other species on the log. Thus, the occurrence of each other predecessor species already on the log became an explanatory variable for the modelled species, and the parameter is interpreted as the species interaction. First, we fit species-level models to experimental data, where competitive interactions were estimated by placing species together and observing their survival and growth (Holmer & Stenlid, 1997; Toljander et al., 2006). Next, we fit models to repeated field inventories where the presence/absence of the species on logs provide colonization-extinction data (Figure 1). Finally, we test if including the species interactions estimated from experimental data as informative priors improve the estimation of species interactions in the field data. Specifically, interactions determined from the experimental data were included as informative priors in the estimation of the interactions in the observational data.
Metacommunity dynamics modelling framework with environmental and species variables
For each successor species, we modeled colonization and extinction dynamics driven by environmental conditions and interactions with other predecessor species (Figure 1) (Kery and Royle, 2020). Each log is a resource unit ‘island’ or independent patch in the meta-community perspective, and the scale at which fungi interact. This clear, well-defined scale helps avoid scale-dependency issues that may hide interactions (Münkemüller et al., 2020).

For all species, models included four environmental variables, both linear and quadratic log diameter and decay stage, and 11 species variables, resulting in 15 colonization variables and 15 extinction variables. The models were fitted using a Bayesian approach, with the prior distributions described below. The estimation was carried out using Markov Chain Monte Carlo simulations with the R library jagsUI (Plummer, 2003; R Core Team, 2018). Three MCMC chains were run for each species, using 50 000 iterations of three chains, of which 10 000 were discarded as burn-in. Model evaluation was carried out through visually examining the posterior chains, and checking if the potential scale reduction statistic (Ȓ) was below 1.1, which indicates chain mixing (Gelman & Hill, 2006). 

Selection of environmental and species interaction variables
Given the large number of environmental and species variables, we used stochastic search variable selection (SSVS) to select important variables for each species model (George & McCulloch, 1993; O’Hara & Sillanpää, 2009). SSVS estimates the probability that individual variables should be included in a regression model. To do so, it uses two normally distributed prior distributions for each variable, but with different variances. More specifically, SSVS uses ‘spike and slab’ prior distributions, and the wide slab distribution (variance σ2 = 1.4, Lemoine, 2019; Northrup and Gerber, 2018) is selected for variables that should be in the model, and the spike distribution (σ2 = 0.01) is selected for the irrelevant variables. If the inclusion of the variable is supported by the data, there will be higher probability for the wide slab prior than the spike prior, and vice versa. We set the prior probability of a wide prior, i.e. a non-zero interspecific interaction, to 0.2 (Mutshinda et al 2009). The MCMC sampling then provides an estimate of the posterior odds of whether a variable should be in the model. Finally, we estimate the amount of evidence for including or not including a variable in the model, the Bayes factor (BF), which is the ratio of the posterior odds and prior odds. We compare the prior odds of a (environmental or species) variable being important (set to 0.2, see above) to the posterior odds estimated by the model. If the quadratic variable was selected by SSVS, we always included the corresponding original linear variable.

George, E.I. & McCulloch, R.E. (1993) Variable Selection via Gibbs Sampling. Journal of the American Statistical Association, 88, 881–889.

Lemoine, N.P. (2019) Moving beyond noninformative priors: why and how to choose weakly informative priors in Bayesian analyses. Oikos, 128, 912–928.

Mutshinda, C.M., O’Hara, R.B. & Woiwod, I.P. (2009) What drives community dynamics? Proceedings of the Royal Society B: Biological Sciences, 276, 2923–2929.

Northrup, J.M. & Gerber, B.D. (2018) A comment on priors for Bayesian occupancy models. PLoS ONE, 13, 1–13.

O’Hara, R.B. & Sillanpää, M.J. (2009) A review of bayesian variable selection methods: What, how and which. Bayesian Analysis, 4, 85–118.


For more information, including model structure, see the article that will soon be published: 

Combining observational and experimental data to estimate environmental and species drivers of fungal metacommunity dynamics
Authors
Hedvig Nenzen1, Helen Moor2, Robert O’Hara3, Mari Jönsson1, Jenni Nordén4, Elisabet Ottosson1, Tord Snäll1
1 SLU Swedish Species Information Centre, Swedish University of Agricultural Sciences, Uppsala, Sweden
2 Swiss Federal Institute for Forest, Snow and Landscape Research WSL, Zurcherstrasse 111, CH-8093 Birmensdorf, Switzerland
3 Department of Mathematical Sciences, Centre for Biodiversity Dynamics, Norwegian University of Science and Technology, N-7491 Trondheim, Norway
4 Norwegian Institute for Nature Research (NINA), Oslo, Norway

Abstract
Understanding the distribution and dynamics of species is central to ecology and important for managing biodiversity. The distributions of many species in metacommunities are determined by many factors, including environmental conditions and interactions between species. Yet, it is difficult to quantify the effect of species interactions on metacommunity dynamics from observational data. We present an approach to estimate the importance of species interactions that combines data from two observational presence-absence inventories (providing colonization-extinction data) with data from species interaction experiments (providing informative prior distributions in the Bayesian framework). We further illustrate the approach on wood-decay fungi that interact within a downed log through competition for resources and space, and facilitate the succession of other species by decaying the wood. Specifically, we estimated the relative importance of species interactions by examining how the presence of a species influenced the colonization and extinction probability of other species. Temporal data on fruit body occurrence of 12 species inventoried twice were jointly analyzed with experimental data from two laboratory experiments aimed to estimate competitive interactions. Both environmental variables and species interactions affect colonization and extinction dynamics. Late-successional fungi had more colonization interactions with predecessor species than early-successional species. We identified several species interactions, and the presence of certain species changed the probability that later-successional species colonized by -81% to 512%. The presence of certain species increased the probability that other species went extinct from a log by 14-61%. Including the informative priors from experimental data added two colonization interactions and one extinction interaction for which the observational field data was inconclusive. However, most species had no detectable interactions, either because they do not interact or because of low species occupancy, meaning data limitation. We show how temporal presence-absence data can be combined with experimental data to identify which species influence the colonization-extinction dynamics of others. Accounting for species interactions in metacommunity models, in addition to environmental drivers, is important because interactions can have cascading effects on other species.
