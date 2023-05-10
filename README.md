# GBIS
The open-source Geodetic Bayesian Inversion Software (GBIS) allows the user to perform the inversion of Interferometric Synthetic Aperture Radar (InSAR) and/or Global Positioning System (GPS) data to estimate deformation source parameters (Bagnardi and Hooper, 2018). The inversion software uses a Markov-chain Monte Carlo algorithm, incorporating the Metropolis-Hastings algorithm [e.g., Hastings, 1970; Mosegaard and Tarantola, 1995], to find the posterior probability distribution of the different source parameters.

GBIS already includes analytical forward models for multiple magmatic sources and fault geometries in elastic half spaces. However, the software architecture allows the user to easily add any other analytical or numerical forward models to calculate displacements at the surface.

GBIS is written in Matlab and uses a series of third-party functions to perform specific steps for data pre-processing, subsampling, and includes forward models from the dMODELS software package [Battaglia et al., 2013].

The software offers a series of pre- and post-processing tools aimed at estimating errors in InSAR data, and at graphically displaying the inversion results.

This version of GBIS has been amended from GBIS_V1.1 (https://comet.nerc.ac.uk/gbis/) and GBIS_V2.1 by members of the Active Tectoncs Group at the University of Leeds.

This is a project of the Center for Observation and Modelling of Volcanoes, Earthquakes and Tectonics (COMET)

References:

Bagnardi M. & Hooper A, (2018). Inversion of surface deformation data for rapid estimates of source parameters and uncertainties: A Bayesian approach. Geochemistry, Geophysics, Geosystems, 19. https://doi.org/10.1029/2018GC007585

Battaglia, M., Cervelli, P. F., & Murray, J. R. (2013). dMODELS: A MATLAB software package for modeling crustal deformation near active faults and volcanic centers. Journal of Volcanology and Geothermal Research, 254, 1-4.

Hastings, W. K. (1970). Monte Carlo sampling methods using Markov chains and their applications. Biometrika, 57(1), 97-109.

Mosegaard, K., & Tarantola, A. (1995). Monte Carlo sampling of solutions to inverse problems. Journal of Geophysical Research: Solid Earth, 100(B7), 12431-12447.
