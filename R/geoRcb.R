.onAttach <- function(lib, pkg) {
  ## Substitute some original functions in geoR
  ## by the ones provided in geoRcb

  ## Function for implanting code into geoR's namespace
  ## code borrowed from utils::assignInNamespace
  replace <- function(x, value, target = 'geoR') {
    ## Namespace of the target
    ## Will be used as the enclosing environment of function x
    ns = getNamespace(target)
    ## Package environment
    ## Where functions bindings are defined
    pe = as.environment(paste('package', target, sep = ':'))

    op <- options()
    if( exists(x, ns, inherits = FALSE) ) {
      if( bindingIsLocked(x, ns) ) {
        unlockBinding(x, ns)
        on.exit({options(warn = -1);
                 lockBinding(x, ns)})
      }
      ## grab the modified function from geoRcb and change its enclosing
      ## environment, as if it was defined within geoR. This is necessary
      ## so that the function can access the full geoR Namespace.
      ## http://adv-r.had.co.nz/Environments.html#function-envs
      value <- get(x, getNamespace('geoRcb'))
      environment(value) <- ns
      ## implant the modified function both in the geoR Namespace
      ## and in the package:geoR environment, which is at this point
      ## already attached
      assign(x, value, envir = ns, inherits = FALSE)

      ## If it is an exported function, its original version was already
      ## attached in the search() path, and we need to substitute it as well
      if( exists(x, pe, inherits = FALSE) ) {
        if( bindingIsLocked(x, pe) ) {
          unlockBinding(x, pe)
          on.exit({options(warn = -1);
                   lockBinding(x, pe)},
                  add = TRUE)
        }
        assign(x, value, envir = pe, inherits = FALSE)
      }
    }
    on.exit(options(op), add = TRUE)

    invisible(NULL)
  }

  ## Functions to be 'implanted' into geoR
  target_funcs <- ls(env = getNamespace('geoRcb'), all.names = TRUE)
  exc.idx <- grep('^\\.__.*?__.*$|^\\.onAttach$|^\\.packageName$',
                  target_funcs)
  target_funcs <- target_funcs[-exc.idx]

  ## Do the job
  for(x in target_funcs) {
    replace(x, x)
  }
}
