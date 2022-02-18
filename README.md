## SWO pseudocode 

The `sims` function is the primary driver, it calls the `pop_est` function and replicates it a user defined number of iterations, it then processes the results out of list form and can save the results as a .csv. 

`pop_est(lfreq = length freq file, cpue = cpue file, samples = NULL, yrs = 2017, strata = NULL)`
            
 - `yrs` are the years to be sampled default >= 2017
 - `samples` is the desired sample size by length  
 - `strata` is a switch to change the comps to year/strata instead of year change to `TRUE` if desired 

1. get complete cases of lengths by species, year, and stratum = `lngs`
    a. there is a `strata` switch to include or exclude
2. filter years from `lfreq` and `cpue` data 
3. calculate length comp
    a. no sampling
        1. comps are calculated by `year`, `species`, `stratum`, `haul`, `sex`, and `length`
    b. reduced sample sizes by sex
        1. filter `lfreq` for only males and females
        2. `uncount` the frequency
        3. sample by `year`, `species`, `stratum`, `haul` at the selected sample size, save as `.new_sexed`
        4. filter out previously sexed fish that are now unsexed, `.new_unsexed`
        5. Join the original unsexed fish from `lfreq` to the `.new_sexed`
            a. *note: this means the overall sample size is reduced, could be computed in another fashion*
        6. comps are calculated by `year`, `species`, `stratum`, `haul`, `sex`, and `length`
4. calculate comp for hauls w/o length samples, `.unk`
5. id hauls w/o lengths, name `.no_length`
6. calculate the population estimate by stratum, `.pop`
7. if there are any hauls w/o lengths joins them to the length comp and compute population size by haul, `.temp`
8. if `strata = TRUE` aggregate abundance by year abd strata, save the output as `new` and the individuals as `removed` in a list
    a. if `strata = NULL` aggregate abundance by year and  output directly to global environment
            

`sims(iters = 1, lfreq = length freq file, cpue = cpue file, strata = NULL, samples = NULL, yrs = 2017, save = NULL)`

 - `iters` is the number of iterations
 - `strata` is a switch to change the comps to year/strata instead of year change to `TRUE` if desired 
 - `samples` is the desired sample size by length  
 - `years` is the min year to sample aka `>=`
 - `save` will save results in the `output` folder, a single file is output if there are no `samples`, otherwise two files are output: a `comps` file and a `removed` file (sexed individuals who are move dto the unsexed population)

1. replicates the `pop_est` function the desired number of iterations


`get_data(data = sim output, id, species = NULL, yrs = NULL)`

 - `id` is the name of the data
 - `species` is a filter 
 - `yrs ` is a filter
 
1. helper function pulls the simulations, splits the results out calculates confidence intervals 

`plot_comp(base_data = og data, sim_data = reduced sample data, species = NULL, yrs = NULL)`

 - species is a filer
 - yrs is a filter
 
1. line plot with 95% ci comparing original population-based comps to reduced sample size comps


`table_comp(base_data = og data, sim_data = reduced sample data, species = NULL, yrs = NULL)`

 - species is a filer
 - yrs is a filter
 
1. creates a table of output values, also output the difference and percent difference by year, or by year & stratum

`plot_comp2(base_data, sim_data, species = NULL, yrs = NULL)`

 - species is a filer
 - yrs is a filter
 
1. density plot to plot the output of `table_comp`


