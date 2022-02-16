## SWO pseudocode 

`pop_est` 

    inputs (lfreq = length freq file
            cpue = cpue file
            samples = NULL
            yrs = 2017
            strata = NULL
            )
            
`yrs` are the years to be sampled default >= 2017  
`samples` is the desired sample size by length  
`strata` is a switch to change the comps to year/strata instead of year change to `TRUE` if desired 

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
            

`sims`

    inputs (iters = 1 
            lfreq = length freq file
            cpue = cpue file
            strata = NULL
            samples = NULL
            yrs = 2017  
            )
`iters` is the number of iterations
`strata` is a switch to change the comps to year/strata instead of year change to `TRUE` if desired 
`samples` is the desired sample size by length  

1. replicate the `pop_est` function the desired number of iterations

`getouts`

    inputs(data = sim output
           type = "comp"
           samples = NULL
           save = NULL
          )

`type` is the result type - only implemented if `!is.null(samples)` pulls either the length comp or the `removed` fish (aka `.new_unsexed`)
`samples` is a flag telling the function to look for a list of results
`save` is a flag to either save the file e.g., `save = "s20.csv"` which will be placed in the `output` folder, if left `NULL` results are placed in the global environment
