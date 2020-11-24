#' Calculate the mean egg counts/rates according to various design types
#'
#' @param simdata The data as simulated by cgpDataSim or clpDataSim
#' @param N A named vector of number of included individuals depending on design type
#' @param count Logical flag to base means on count data or the underlying rates
#' @param log_constant A constant to add to the count data before calculating geometric means (ignored if count==FALSE)
#' @param screen_threshold The threshold count on which to screen individuals
#'
#' @return A data frame of summary statistics
#' @examples
#' data <- cgpDataSim(10^3, 600, 0.1, 100, 1, 1, 1, 0, true_prevalence=0.8)
#' means <- design_means(data, second_slide_cost = c(0.1, 0.621, 1))
#' library("dplyr")
#' means %>%
#'   group_by(Design, SecondSlideCost) %>%
#'   summarise(N = mean(N), Budget = mean(ScreenBudget + SampleBudget),
#'     Bias = median(ArithmeticEfficacy - TrueArithmetic),
#'     SD = sd(ArithmeticEfficacy - TrueArithmetic))
#'
#' @importFrom tidyr replace_na expand_grid
#' @import dplyr
#'
#' @export
design_means <- function(simdata, design = c('NS','NS2','NS3','SS','SSR1','SSR2','SSR3'), budget=600, second_slide_cost = 0.621, max_screen = 0.9, count=TRUE, log_constant=if(count) 1 else 0, screen_threshold = 0)
#design_means <- function(simdata, N = c(NS=600, SS=109, SSR1=100, SSR2=92, SSR3=95), count=TRUE, log_constant=if(count) 1 else 0, screen_threshold = 0)
{
  # TODO: check simdata is valid etc

  stopifnot(all(design %in% c('NS','NS1','NS2','NS3','SS','SSR1','SSR2','SSR3')))

  # Allow second_slide_cost to be vectorised where relevant:
  designcombo <- expand_grid(design = design, budget=budget, second_slide_cost = NA_real_) %>%
    filter(! design %in% c("NS3","SSR3")) %>%
    bind_rows( expand_grid(design = design[design %in% c("NS3","SSR3")], budget=budget, second_slide_cost = second_slide_cost))

  res <- bind_rows(lapply(seq_len(nrow(designcombo)), function(i) getmeans_slideN(simdata, design=designcombo$design[i], budget=designcombo$budget[i], second_slide_cost = designcombo$second_slide_cost[i], max_screen = max_screen, count=count, log_constant=log_constant, screen_threshold=screen_threshold)))

  return(res)
}



# Helper function:
getmeans_slideN <- function(simdata, design='NS', budget=600, second_slide_cost = 0.621, max_screen = 0.9, count=TRUE, log_constant=if(count) 1 else 0, screen_threshold = 0)
{

  # TODO: ensure that simdata$CVbetween:CVreduction doesn't vary

  # Select either count or lambda and add the constant:
  if(count){
    simdata <- simdata %>%
      mutate(ScreenUsing = ScreenCount, PreUsing = PreCount, PostUsing1a = PostCount1a, PostUsing1b = PostCount1b, PostUsing2 = PostCount2)
  }else{
    simdata <- simdata %>%
      mutate(ScreenUsing = ScreenSlide, PreUsing = PreSlide, PostUsing1a = PostSlide1a, PostUsing1b = PostSlide1b, PostUsing2 = PostSlide2)
  }

  stopifnot(is.numeric(max_screen) && length(max_screen)==1 && max_screen >0 && max_screen <1)

  # Then add columns for the best estimated geometric and arithmetic mean reductions based on the NS type:
  simdata <- simdata %>%
    group_by(Reduction) %>%
    mutate(BestArithmetic = 100*(1-mean(PostSlide1a)/mean(PreSlide))) %>%
    mutate(BestGeometric = 100*(1-exp(mean(log(PostSlide1a))-mean(log(PreSlide))))) %>%
    ungroup()


  # Ensure that the budget is equal between communities:
  community_budget <- budget / length(unique(simdata$Community))

  # Ensure a sufficient amount of data is available:
  ns <- simdata %>% count(Replicate, Reduction, Community)
  stopifnot(all(ns$n >= community_budget))


  # Helper function to round up or down randomly (not being used currently):
  # ceiloor <- function(x) floor(x) + rbinom(length(x),1,x%%1)

  # Then subselect data according to design:
  lc <- log_constant
  if(design%in%c('NS','NS1')){
    # We simply take the first community_budget/2 individuals per simulation:
    subsampled <- simdata %>%
      group_by(Replicate, Reduction, Community) %>%
      mutate(ScreenBudget = 0, SampleBudget = (1:n())*2) %>%
      filter((ScreenBudget + SampleBudget) <= community_budget) %>%
      ungroup() %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = NA)
  }else if(design=='NS2'){
    # We simply take the first community_budget/3 individuals per simulation:
    subsampled <- simdata %>%
      group_by(Replicate, Reduction, Community) %>%
      mutate(ScreenBudget = 0, SampleBudget = (1:n())*3) %>%
      filter((ScreenBudget + SampleBudget) <= community_budget) %>%
      ungroup() %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = PostUsing2)
  }else if(design=='NS3'){
    # We simply take the first community_budget/(2+ssc) individuals per simulation:
    subsampled <- simdata %>%
      group_by(Replicate, Reduction, Community) %>%
      mutate(ScreenBudget = 0, SampleBudget = (1:n())*(2+second_slide_cost)) %>%
      filter((ScreenBudget + SampleBudget) <= community_budget) %>%
      ungroup() %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = PostUsing1b)
  }else if(design=='SS'){
    # We do new pre samples until we have exhausted the budget:
    subsampled <- simdata %>%
      group_by(Replicate, Reduction, Community) %>%
      mutate(ScreenBudget = cumsum(PreUsing <= screen_threshold), SampleBudget = 2*cumsum(PreUsing > screen_threshold)) %>%
      filter((ScreenBudget + SampleBudget) <= community_budget) %>%
      filter(PreUsing > screen_threshold) %>%
      ungroup() %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = NA)
  }else if(design=='SSR1'){
    # We do new screening samples until we have exhausted the budget:
    subsampled <- simdata %>%
      group_by(Replicate, Reduction, Community) %>%
      mutate(ScreenBudget = 1:n(), SampleBudget = 2*cumsum(ScreenUsing > screen_threshold)) %>%
      filter((ScreenBudget + SampleBudget) <= community_budget) %>%
      filter(ScreenUsing > screen_threshold) %>%
      ungroup() %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = NA)
  }else if(design=='SSR2'){
    # We do new screening samples until we have exhausted the budget:
    subsampled <- simdata %>%
      group_by(Replicate, Reduction, Community) %>%
      mutate(ScreenBudget = 1:n(), SampleBudget = 3*cumsum(ScreenUsing > screen_threshold)) %>%
      filter((ScreenBudget + SampleBudget) <= community_budget) %>%
      filter(ScreenUsing > screen_threshold) %>%
      ungroup() %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = PostUsing2)
  }else if(design=='SSR3'){
    # We do new screening samples until we have exhausted the budget:
    subsampled <- simdata %>%
      group_by(Replicate, Reduction, Community) %>%
      mutate(ScreenBudget = 1:n(), SampleBudget = (2+second_slide_cost)*cumsum(ScreenUsing > screen_threshold)) %>%
      filter((ScreenBudget + SampleBudget) <= community_budget) %>%
      filter(ScreenUsing > screen_threshold) %>%
      ungroup() %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = PostUsing1b)
  }else{
    stop('Unrecognised design strategy')
  }

  ## If all pre-tx samples are zero and/or we have used more than screen_max of the budget on screening then throw out the community:
  budgets <- subsampled %>%
    group_by(Replicate, TruePrev, ObsPrev, Reduction, Community) %>%
    summarise(PreMean = mean(Pre), ScreenBudget = max(ScreenBudget), SampleBudget = max(SampleBudget), ScreenProp = ScreenBudget / (ScreenBudget+SampleBudget), .groups='drop') %>%
    filter(ScreenProp <= max_screen, PreMean > 0) %>%
    group_by(Replicate, TruePrev, ObsPrev, Reduction) %>%
    mutate(Communities = length(unique(Community))) %>%
    ungroup()

  ## Then calculate summary statistics:
  sumstats <- subsampled %>%
    select(-PreMean, -ScreenBudget, -SampleBudget, -TruePrev, -ObsPrev) %>%
    inner_join(budgets, by = c("Replicate", "Community", "Reduction")) %>%
    group_by(Replicate, TruePrev, ObsPrev, Reduction, TrueGeometric, TrueArithmetic, BestArithmetic, BestGeometric, OverallMean, ScreenBudget, SampleBudget, ScreenProp, Communities) %>%
    summarise(N = n(), Count = count, ArithmeticEfficacy = 100*(1-mean(c(Post1, Post2), na.rm=TRUE)/mean(Pre)), GeometricEfficacy = 100*(1-exp(mean(log(c(Post1, Post2)+lc), na.rm=TRUE)-mean(log(Pre+lc)))), .groups='drop')

  ## If we have completely removed a replicate (all communities gone) then re-add it:
  sumstats <- subsampled %>%
    group_by(Replicate, Reduction, TrueGeometric, TrueArithmetic, BestArithmetic, BestGeometric, OverallMean) %>%
    summarise(.groups='drop') %>%
    full_join(sumstats, by=c("Replicate", "Reduction", "TrueGeometric", "TrueArithmetic", "BestArithmetic", "BestGeometric", "OverallMean")) %>%
    mutate(Design = design, Budget = budget, SecondSlideCost = second_slide_cost,
           Communities = replace_na(Communities, 0))

  if(any(sumstats$Communities == 0)) warning("One or more replicate found nobody to sample")

  return(sumstats)
}



##### NOT USING:
# Helper function:
getmeans_fixedN <- function(simdata, N, design='NS', count=TRUE, log_constant=if(count) 1 else 0, screen_threshold = 0)
{

  # TODO: ensure that cv, mean and reduction doesn't vary

  # Select either count or lambda and add the constant:
  if(count){
    simdata <- simdata %>%
      mutate(ScreenUsing = ScreenCount, PreUsing = PreCount, PostUsing1a = PostCount1a, PostUsing1b = PostCount1b, PostUsing2 = PostCount2)
  }else{
    simdata <- simdata %>%
      mutate(ScreenUsing = ScreenSlide, PreUsing = PreSlide, PostUsing1a = PostSlide1a, PostUsing1b = PostSlide1b, PostUsing2 = PostSlide2)
  }

  # Then add columns for the best estimated geometric and arithmetic mean reductions based on the NS type:
  simdata <- simdata %>%
    group_by(Reduction) %>%
    mutate(Prevalence = sum(ScreenCount > screen_threshold) / n()) %>%
    mutate(BestArithmetic = 100*(1-mean(PostSlide1a)/mean(PreSlide))) %>%
    mutate(BestGeometric = 100*(1-exp(mean(log(PostSlide1a))-mean(log(PreSlide))))) %>%
    ungroup()

  # Then subselect data according to design:
  lc <- log_constant
  if(design=='NS'){
    subsampled <- simdata %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = NA)
  }else if(design=='SS'){
    subsampled <- simdata %>%
      filter(PreCount > screen_threshold) %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = NA)
  }else if(design=='SSR1'){
    subsampled <- simdata %>%
      filter(ScreenCount > screen_threshold) %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = NA)
  }else if(design=='SSR2'){
    subsampled <- simdata %>%
      filter(ScreenCount > screen_threshold) %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = PostUsing2)
  }else if(design=='SSR3'){
    subsampled <- simdata %>%
      filter(ScreenCount > screen_threshold) %>%
      mutate(Pre = PreUsing, Post1 = PostUsing1a, Post2 = PostUsing1b)
  }else{
    stop('Unrecognised design strategy')
  }

  # Ensure a balance between communities:
  C <- length(unique(simdata$Community))
  NC <- ceiling(N/C)

  # Then calculate summary statistics:
  sumstats <- subsampled %>%
    group_by(Replicate, Reduction, Community) %>%
    slice(1:min(n(), NC)) %>%
    #sample_n(min(NC, n())) %>%
    ungroup() %>%
    group_by(Replicate, Prevalence, Reduction, TrueGeometric, TrueArithmetic, BestArithmetic, BestGeometric, OverallMean) %>%
    summarise(Design = design, N = n(), Count = count, ArithmeticEfficacy = 100*(1-mean(c(Post1, Post2), na.rm=TRUE)/mean(Pre)), GeometricEfficacy = 100*(1-exp(mean(log(c(Post1, Post2)+lc), na.rm=TRUE)-mean(log(Pre+lc))))) %>%
    ungroup()

  if(any(sumstats$N < N))
    warning('One or more maxN was insufficient to obtain the specified N')

  return(sumstats)
}
##### \NOT USING
