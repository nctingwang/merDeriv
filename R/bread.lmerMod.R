bread.lmerMod <- function (x, ...){  
  ## get all elements by getME and exclude multiple random effect models.
  if (!is(x, "lmerMod")) stop("bread.lmerMod() only works for lmer() models.")
  
  parts <- getME(x, "ALL")
  if (length(parts$l_i) > 1) stop("Multiple cluster variables detected. Robust SEs are unavailable.")
  
  dotdotdot <- list(...)
  if("full" %in% names(dotdotdot) | "information" %in% names(dotdotdot)){
    full <- dotdotdot$full
    information <- dotdotdot$information
  } else {
    full <- FALSE
    information <- "expected"
  }
  if(!('full' %in% c("TRUE", "FALSE"))) stop("invalid 'full' argument supplied")
  if(!('information' %in% c("expected", "observed"))) stop("invalid 'information' argument supplied")
  object <- x
  return(vcov.lmerMod(object, full = full, information = information) * parts$l_i)
}
