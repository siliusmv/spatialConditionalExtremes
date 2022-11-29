
#' @export
data_dir = function() file.path(here::here(), "inst", "extdata")
#' @export
downloads_dir = function() file.path(data_dir(), "downloads")
#' @export
image_dir = function() file.path(data_dir(), "images")
#' @export
results_dir = function() file.path(data_dir(), "results")
#' @export
cgeneric_dir = function() file.path(here::here(), "cgeneric")

#' Compute the logarithm of the sum of the exponents of a vector x.
#' This function is usable in settings where components x_i of x are so large/big that
#' exp(x_i) = Inf, or exp(x_i) = 0
#' @export
log_sum_exp = function(x, na.rm = FALSE) {
  m = max(x, na.rm = na.rm)
  if (is.infinite(m) && m < 0) {
    res = -Inf
  } else {
    res = m + log(sum(exp(x - m), na.rm = na.rm))
  }
  res
}

#' Compute the logarithm of the mean of the exponents of a vector x.
#' This function is usable in settings where components x_i of x are so large/big that
#' exp(x_i) = Inf, or exp(x_i) = 0
#' @export
log_mean_exp = function(x, na.rm = FALSE) {
  log_sum_exp(x, na.rm) - log(sum(!is.na(x)))
}

#' Compute the Euclidean distance between x and y, where
#' x is either a row vector, or a matrix of row vectors,
#' and y is a matrix of row vectors
#' @export
dist_euclid = function(x, y) {
  if (!is.matrix(y)) y = matrix(y, ncol = 1)
  if (is.matrix(x)) {
    res = dist_euclid_mat_mat(x, y)
  } else {
    res = dist_euclid_vec_mat(x, y)
  }
  res
}

#' Execute a shell script, and provide feedback if it crashes
#' @export
execute_shell_script = function(command, args = character(), ...) {
  output = system2(command, args, ...)
  success = (output == 0)
  if (!success) {
    formatted_args = paste(args, collapse = " ")
    stop("shell script ", command,
         " gave error code ", output,
         " with arguments ", formatted_args)
  }
  0
}

#' Call a makefile used for compiling and linking cgeneric scripts,
#' in order to use them with R-INLA
#' @export
make_cgeneric = function(cmd) {
  current_path = getwd()
  on.exit(setwd(current_path))
  setwd(cgeneric_dir())
  execute_shell_script("make", cmd)
}

#' Create a progress bar for tracking long-running processes.
#' This is a thin wrapper around the progress package
#' @export
progress_bar = function(n) {
  pb = progress::progress_bar$new(
    format = ":percent [:bar] time elapsed: :elapsedfull, eta: :eta",
    total = n, width = 70, clear = FALSE)
  # pb$tick() throws an error if we tick too many times, which can potentially stop
  # a script that is otherwise working fine. We don't want that to happen,
  # so we change the tick function slightly
  res = list(
    terminate = pb$terminate,
    tick = function(...) tryCatch(pb$tick(...), error = function(e) NULL))
  res
}

#' Turn an R plot into a beautiful pdf made by LaTeX and TikZ,
#' using the tikzDevice package
#' @export
tikz_plot = function(file, plot = NULL, expression = NULL, ...) {
  # Ensure that you are on an operating system that you have tested
  operating_system = Sys.info()[["sysname"]]
  if (operating_system == "Windows") {
    proceed = readline(paste("This function was written on a Mac,",
                             "I have no idea if it will work on Windows.",
                             "Proceed? (y/n) "))
    if (proceed != "y") return()
  }

  # Create a temporary file for the tikz-output
  tmp = tempfile(tmpdir = getwd())
  # Clean up after yourself on early interrupt
  on.exit(suppressWarnings(file.remove(tmp)), add = TRUE)

  # Extract default tex usepackages and add the bm package for bold greek letters
  opt = options()
  on.exit(options(opt)) #Reset global options on exit
  tikzDevice::setTikzDefaults(overwrite = FALSE)
  tex_packages = options()$tikzLatexPackages
  if (!any(grepl("usepackage\\{bm\\}", tex_packages))) {
    tex_packages = c(tex_packages, "\\usepackage{bm}\n")
  }
  if (!any(grepl("usepackage\\{amsmath\\}", tex_packages))) {
    tex_packages = c(tex_packages, "\\usepackage{amsmath}\n")
  }

  # Open a device for creating a tex-file
  tikzDevice::tikz(tmp, standAlone = TRUE, packages = tex_packages, ...)
  # Call dev.off() on exit in case of interruptions
  current_device = dev.cur()
  on.exit(dev.off(current_device))

  # Plot something into the tex-file
  if (!is.null(plot)) {
    if (any(class(plot) %in% c("gg", "ggplot", "patchwork"))) {
      print(plot)
    } else {
      for (p in plot) print(p)
    }
  } else {
    eval(substitute(expression), envir = parent.frame())
  }

  # Finish the creation of the tex-file
  dev.off()

  # Compile to pdf using lualatex
  system2("lualatex", shQuote(tmp))

  # Copy pdf file to final destination
  file.copy(paste0(tmp, ".pdf"), file, overwrite = TRUE)

  # Clean up all temporary files
  tmp_filename = tail(strsplit(tmp, "/")[[1]], 1)
  files_to_clean = grep(tmp_filename, list.files(full.names = TRUE), value = TRUE)
  file.remove(files_to_clean)
}

dist_euclid_vec_mat = function(x, y) {
  stopifnot(length(x) == ncol(y))
  apply(y, 1, function(y) sqrt(sum((x - y)^2)))
}

dist_euclid_mat_mat = function(x, y) {
  apply(x, 1, dist_euclid_vec_mat, y = y)
}
