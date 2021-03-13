# eggSim
eggSim is an R package developed by Matt Denwood and Luc Coffeng to assess different survey designs for faecal egg count reduction tests, in terms of their cost efficiency, precision and expected statistical bias. The code is open-source and freely available for use, although it is currently under active development so expect new features to be added over time. We also welcome feedback and contributions to the project.

## Website

For a background to our work along with simple installation and usage instructions see our website:  [https://www.fecrt.com/surveys](https://www.fecrt.com/surveys)

## TODO

- Fix the "object 'second_slide_cost' not found" error when running in parallel under Windows (not exported to the PSOCK cluster?)
- Re-write the data simulation/analysis in C++
- Harmonise the lognormal data simulation with the gamma data simulation
- autoplot and summary methods
- Expand the vignettes
