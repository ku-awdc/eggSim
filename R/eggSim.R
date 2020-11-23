#' Simulate and summarise egg reduction rate data
#'
#' @param R The number of replicate datasets
#' @param N The maximum number of individuals in total
#' @param reduction The arithmetic mean reduction (may be vectorised)
#' @param community_mean A vector of arithmetic mean pre-treatment EPG in each community
#' @param cv_between A vector of CV reflecting variation in EPG between individuals in each community
#' @param cv_within Day-to-day variation in EPG within an individual
#' @param cv_slide Variation between faecal samples from the same individual and day
#' @param cv_reduction Variation in efficacy between individuals
#' @param edt The egg detection threshold (24 EPG is standard for Kato-Katz)
#' @param N A named vector of number of included individuals depending on design type
#' @param count Logical flag to base means on count data or the underlying rates
#' @param log_constant A constant to add to the count data before calculating geometric means (ignored if count==FALSE)
#' @param screen_threshold The threshold count on which to screen individuals
#' @param maxN The maximum number of simulated individuals (on which subsampling occurs)
#'
#' @return A data frame containing the simulated data
#'
#' @examples
#' means <- eggSim(c(0.2,0.1,0.05), cv_reduction=0)
#'
#' @importFrom parallel detectCores makeForkCluster makePSOCKcluster clusterSetRNGStream parLapply stopCluster
#' @importFrom pbapply pblapply
#' @importFrom purrr walk
#'
#' @export
eggSim <- function(reduction, budget=600, second_slide_cost = 0.621, max_screen = 0.9, community_mean=c(24, 48), cv_between=c(1.5), cv_within=0.75, cv_slide=0.25, cv_reduction=0, count=TRUE, log_constant=if(count) 1 else 0, screen_threshold = 0, grams=1/24, R=10^3, summarise = TRUE, type="gamma", parallelise=TRUE, max_vec = 5e6)
{

  st <- Sys.time()

  ## TODO: allow vectorisation over community_mean, reduction and budget
  ## community_mean on separate jobs (not vectorised internally)
  ## reduction split over cores (vectorised internally)
  ## budget (can be done from the same simulated data)

  if(!is.matrix(community_mean)){
    ## Shortcut to allow vectorising mean of a single community:
    if(length(cv_between) == 1){
      dim(community_mean) <- c(length(community_mean), 1)
    }else{
    ## Shortcut in case mean is vectorised by community only:
      stopifnot(length(community_mean) != length(cv_between))
      dim(community_mean) <- c(1, length(community_mean))
    }
  }
  stopifnot(ncol(community_mean) == length(cv_between))
  meanindex <- seq_along(community_mean)
  combos <- expand_grid(reduction, meanindex)

  # Determine the minimum number of splits so that all means
  # are on different jobs, and so that no individual job creates
  # a data frame with more than max_vec rows
  rows_per_red <- R*max(budget)*ncol(community_mean)
  if(rows_per_red > max_vec){
    max_vec <- rows_per_red
    warning("The specified max_vec was insufficient: increasing to ", max_vec)
  }
  preds <- ((rows_per_red*length(reduction)) %/% max_vec) +1
  minsplits <- length(meanindex) * preds
  ## TODO: there is a problem with combining reductions
  ## BUT it is faster to run serially (with less RAM) anyway
  minsplits <- length(meanindex) * length(reduction)

  ## TODO: vectorise budget in design_means
  stopifnot(length(budget)==1)

  ## TODO: bring lognormal data option back!
  stopifnot(length(type) == 1 && type %in% "gamma")

  # TODO: complete argument checks
  stopifnot(length(cv_within) == 1 && cv_within >= 0.0)
  stopifnot(length(cv_slide) == 1 && cv_slide >= 0.0)
  stopifnot(length(cv_reduction) == 1 && cv_reduction >= 0.0)
  stopifnot(length(grams) == 1 && grams >= 0.0)

  cl <- NULL
  cores <- getOption("mc.cores", detectCores())

  if(inherits(parallelise, "cluster")){
    cl <- parallelise
    cores <- length(as.character(cl))
    parallelise <- TRUE
  }

  if(is.numeric(parallelise)){
    stopifnot(length(parallelise)==1, parallelise > 0, parallelise%%1 == 0)
    cores <- parallelise
    parallelise <- TRUE
  }

  stopifnot(is.logical(parallelise), length(parallelise)==1)
  if(parallelise && length(reduction) > 1){
    combos <- split(combos, 1+ (seq_len(nrow(combos))-1)%%max(minsplits, cores))

    if(!inherits(cl, "cluster")){
      cores <- min(cores, length(combos))
      if(.Platform$OS.type=="unix"){
        cl <- makeForkCluster(cores)
      }else{
        cl <- makePSOCKcluster(cores)
      }
      clusterSetRNGStream(cl, NULL)
      on.exit(stopCluster(cl))
    }
    stopifnot(inherits(cl, "cluster"))

  }else{
    combos <- combos <- split(combos, 1+ (seq_len(nrow(combos))-1)%%minsplits)
  }

  # Ensure that no job has more than 1 unique meanindex:
  combos %>%
    walk(~ if(any(.x$meanindex != .x$meanindex[1])) stop("Logic error in splitting jobs"))

  # Worker function:
  fun <- function(comb, silent=TRUE){

    commu <- as.numeric(community_mean[comb$meanindex,,drop=TRUE])
    red <- comb$reduction

    # First simulate gamma and lognormal data:
    if(!silent) cat('Simulating gamma data...\n')
    simdatagp <- cgpDataSim(R, max(budget), red, commu, cv_between, cv_within, cv_slide, cv_reduction, grams=grams)
    #if(!silent) cat('Simulating lognormal data...\n')
    #simdatalp <- clpDataSim(R, max(budget), red, commu, cv_between, cv_within, cv_slide, cv_reduction, grams=grams)

    # Then summarise:
    if(!silent) cat('Summarising gamma data...\n')
    meansgp <- design_means(simdatagp, budget=budget, second_slide_cost=second_slide_cost, max_screen=max_screen, count=count, log_constant=log_constant, screen_threshold=screen_threshold) %>% mutate(ComparisonArithmetic = TrueArithmetic, ComparisonGeometric = BestGeometric, IncludeNS = "Arithmetic", Data = if(count) 'Gamma-Poisson' else 'Gamma')

    #if(!silent) cat('Summarising lognormal data...\n')
    #meanslp <- design_means(simdatalp, N=N, count=count) %>% mutate(ComparisonArithmetic = BestArithmetic, ComparisonGeometric = TrueGeometric, IncludeNS = "Geometric", Data = if(count) 'Lognormal-Poisson' else 'Lognormal')

    #return(bind_rows(meansgp,meanslp))
    return(meansgp)

  }

  if(length(combos)==1){
    cat("Running a single set of simulations...\n")
    output <- fun(combos[[1]], silent=FALSE)
  }else if(!parallelise){
    cat("Running ", length(combos), " sets of simulations in serial...\n", sep="")
    output <- pblapply(combos, fun, silent=TRUE, cl=NULL) %>%
      bind_rows()
  }else{
    cat("Running ", length(combos), " sets of simulations over ", min(cores, length(combos)), " parallel cores...\n", sep="")
    output <- pblapply(combos, fun, silent=TRUE, cl=cl) %>%
      bind_rows()
  }

  if(!summarise) return(output)


  # Then join data frames:
  cat('Summarising output...\n')
  means <- output %>%
    select(Design, Prevalence, ComparisonArithmetic, ComparisonGeometric, IncludeNS, Data, OverallMean, Count, ScreenProp, N, Communities, ArithmeticEfficacy, GeometricEfficacy) %>%
    gather(Type, Efficacy, -Design, -ComparisonArithmetic, -ComparisonGeometric, -IncludeNS, -Data, -OverallMean, -Count, -Communities, -ScreenProp, -N, -Prevalence) %>%
    mutate(Type = gsub('Efficacy','',Type), Target = ifelse(Type=='Arithmetic', ComparisonArithmetic, ComparisonGeometric)) %>%
    mutate(Set = paste0(Data, ' - ', Type), Cheating = Design=='NS' & IncludeNS!=Type) %>%
    group_by(Design, Set, Data, Type, Target, OverallMean) %>%
    mutate(Success = sum(Communities > 0) / n()) %>%
    ungroup() %>%
    filter(Communities > 0) %>%
    group_by(Design, Set, Data, Type, Cheating, OverallMean, Target, Success) %>%
    summarise(Prevalence = mean(Prevalence), ScreenProp = mean(ScreenProp), MeanN = mean(N), Bias = mean(Efficacy - Target), MedianBias = median(Efficacy - Target), MeanRatio = mean(Efficacy/Target), ReductionRatio = mean((1-Efficacy/100) / (1-Target/100)), LCI = Bias - 1.96*sd(Efficacy - Target)/sqrt(n()), UCI = Bias + 1.96*sd(Efficacy - Target)/sqrt(n()), Variance = var(Efficacy), VarianceMeanRatio = var(Efficacy/Target)) %>%
    ungroup()

  # Remove bias estimates where it is cheating:
  means$Bias[means$Cheating] <- NA
  means$LCI[means$Cheating] <- NA
  means$UCI[means$Cheating] <- NA
  means$Cheating <- NULL

  # Filter out irrelevant data types:
  if(type=="gamma"){
    means <- means %>%
      filter(Type=="Arithmetic")
  }

  cat("Done (time elapsed: ", round(as.numeric(Sys.time()-st, units="mins"),1), " minutes)\n", sep="")

  return(means)

}

