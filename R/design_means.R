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
#'
#' @export
design_means <- function(simdata, N = c(NS=600, SS=109, SSR1=100, SSR2=92, SSR3=95), count=TRUE, log_constant=if(count) 1 else 0, screen_threshold = 0)
{
  # TODO: check simdata is valid etc

  stopifnot(all(names(N) %in% c('NS','SS','SSR1','SSR2','SSR3')))
  res <- bind_rows(lapply(1:length(N), function(i) getmeans(simdata, N[i], names(N)[i], count, log_constant, screen_threshold)))

  return(res)
}

# Helper function:
getmeans <- function(simdata, N, design='NS', count=TRUE, log_constant=if(count) 1 else 0, screen_threshold = 0)
{

  # TODO: ensure that cv stuff doesn't vary

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
    #slice(1:NC) %>%
    sample_n(min(NC, n())) %>%
    ungroup() %>%
    group_by(Replicate, Reduction, TrueGeometric, TrueArithmetic, BestArithmetic, BestGeometric, OverallMean) %>%
    summarise(Design = design, N = n(), Count = count, ArithmeticEfficacy = 100*(1-mean(c(Post1, Post2), na.rm=TRUE)/mean(Pre)), GeometricEfficacy = 100*(1-exp(mean(log(c(Post1, Post2)+lc), na.rm=TRUE)-mean(log(Pre+lc))))) %>%
    ungroup()

  if(any(sumstats$N < N))
    warning('One or more maxN was insufficient to obtain the specified N')

  return(sumstats)
}

