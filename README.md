# South African TB Hotspot Model

## Synopsis

This R project contains a matematical model to study the effect of certain public health interventions on tuberculosis (TB) dynamics between the general community and a presumed hotpot. The current implementation of this model is deterministic. 

## Motivation

This project was started as an assignment in the 2015 offering of Infectious Disease Dynamics at The Johns Hopkins University Bloomberg School of Public Health. The course instructors were Justin Lessler, PhD and Derek Cummings, PhD. This repository seeks to expand the work created and submitted under the requirements for the assignment. 

The model is based on the TB model described by Dowdy et al. in "Heterogeneity in tuberculosis transmission and the role of geographic hotspots in propagating epidemics."

## Description of the System
![](https://raw.githubusercontent.com/alvinthomas/SA.TBHotspot/master/system_diagram.png)

This model defines two populations: the host (or community) and the hotspot. Within each population, there people can exist in one of 10 compartments. Half of these compartments are for individuals with HIV and the other half are for individuals without HIV. The five types of compartments are susceptibles, latent (recently infected), latent (remotely infected), active infection, and recovered. 

The current model uses a steady-state design. That is to say, the number of deaths always equal the number of births. The system can be defined by the following system of differential equations.

![](https://raw.githubusercontent.com/alvinthomas/SA.TBHotspot/master/system_of_eqns.png)

Additionally, I use the following system set of parameter equations. 
![](https://raw.githubusercontent.com/alvinthomas/SA.TBHotspot/master/system_of_eqns2.png)

## Installation

TODO

## License
Copyright (c) 2015 Alvin Thomas

Please refer to the LICENSE file containing the MIT License (MIT) for this project.

