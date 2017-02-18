bread.lmerMod <- function (object, full = TRUE, ...){  
  ## get all elements by getME and exclude multiple random effect models.
  if (!is(object, "lmerMod")) stop("bread.lmerMod() only works for lmer() models.")
  
  parts <- getME(object, "ALL")
  if (length(parts$l_i) > 1) stop("Multiple cluster variables detected. Robust SEs are unavailable.")
  
  return(vcov.lmerMod(object, full = full) * parts$l_i)
}
