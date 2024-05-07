# Patter.jl

**Particle filters, smoothers and samplers for animal movement modelling in [`Julia`](https://julialang.org)**

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Documenter](https://github.com/edwardlavender/Patter.jl/actions/workflows/Documenter.yml/badge.svg)](https://github.com/edwardlavender/Patter.jl/actions/workflows/Documenter.yml)
[![Runtests](https://github.com/edwardlavender/Patter.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/edwardlavender/Patter.jl/actions/workflows/Runtests.yml)
[![Coverage](https://codecov.io/gh/edwardlavender/Patter.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/edwardlavender/Patter.jl)

[`Patter.jl`](https://github.com/edwardlavender/Patter.jl) is a [`Julia`](https://julialang.org) package that provides particle filters, smoothers and samplers for animal movement modelling, with a focus on passive acoustic telemetry systems. [`Patter.jl`](https://github.com/edwardlavender/Patter.jl) is heavily based on the [`ParticleFish`](https://github.com/scheidan/ParticleFish.jl) package developed by [Andreas Scheidegger](https://www.eawag.ch/de/ueber-uns/portraet/organisation/mitarbeitende/profile/andreas-scheidegger/show/). The package forms the backend for the [`patter`](https://github.com/edwardlavender/patter) [`R`](https://www.r-project.org) package.

# Features

[`Patter.jl`](https://github.com/edwardlavender/Patter.jl) is designed to reconstruct fine-scale movement paths and emergent patterns of space use for tagged animals. Powerful, Monte-Carlo algorithms known as particle filters, smoothers and samplers are used for this purpose. The essential functions are `particle_filter()`, `particle_twofilter()` and `particle_sampler()`. `particle_filter()` is the particle filter. This simulates the possible locations of an individual moving forwards in time, accounting for all of the data (such as detections at receivers and depth measurements) up to each time point and the animal’s movement. `particle_twofilter()` and `particle_sampler()` refine outputs from the particle filter, accounting for all of the data _up to and after_ each time point. 

# Installation

* Install [`Julia`](https://julialang.org) ≥ v.1.9;

* Install [`Patter.jl`](https://github.com/edwardlavender/Patter.jl):
    - Use `]` in the Julia REPL to open the package manager;
    - Use `add https://github.com/edwardlavender/Patter.jl` to install [`Patter.jl`](https://github.com/edwardlavender/Patter.jl);
    - Alternatively, use 

You can also [`Patter.jl`](https://github.com/edwardlavender/Patter.jl) via the [`patter`](https://github.com/edwardlavender/patter) [`R`](https://www.r-project.org) wrapper.

# Functions

## Simulation 

To simulate datasets, use:

* `simulate_path_walk()` to simulate a movement path from a walk model;
* `simulate_observations()` to simulate observations, such as detections at receivers;

## Real-world datasets

To collate real-world datasets for real-world analyses, use:

* `assemble_yobs()` to assemble a hash-table;

This function expects a `Vector` of `DataFrame`s, one for each data type, that comprise a timeline of observations and associated model parameters, and a corresponding `Vector` of observation model (`ModelObs`) subtypes. 

## Modelling

To simulate initial states (i.e., locations) for the particle filter, use: 
* `simulate_states_init()` to simulate states uniformally across an area;

To define a movement model, see:
* `ModelMove` to create a movement model instance;

For available observation models, see:
* `ModelObs`

To implement the particle filter, use:
* `particle_filter()` to run the filter;

# Usage 

TO DO

# Citation

Lavender, E. et al. (in prep). Particle filters, smoothers and samplers for animal movement modelling.

--- 

Please note that [`Patter.jl`](https://github.com/edwardlavender/Patter.jl) is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.