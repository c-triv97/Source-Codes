# ---- Write ---------------------------------------------------------------
write_requirements <- function(file = "requirements.txt") {
  pkgs <- sessionInfo()$otherPkgs          # attached, non-base packages
  if (length(pkgs) == 0) {
    stop("No non-base packages are attached.", call. = FALSE)
  }

  versions <- vapply(pkgs, function(x) x$Version, character(1))
  lines <- paste0(names(pkgs), "==", versions)

  writeLines(sort(lines), file)
  message("Wrote ", length(lines), " packages to ", file)
  invisible(lines)
}

# ---- Read + load ---------------------------------------------------------
load_requirements <- function(
  file = "requirements.txt",
  install_missing = TRUE,
  check_version = TRUE
) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)

  lines <- trimws(readLines(file, warn = FALSE))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]

  parts <- strsplit(lines, "==", fixed = TRUE)
  pkgs  <- vapply(parts, `[`, character(1), 1)
  vers  <- vapply(parts, function(x) if (length(x) > 1) x[2] else NA_character_,
                  character(1))

  # Install anything missing (BiocManager handles both CRAN and Bioconductor)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0 && install_missing) {
    message("Installing: ", paste(missing, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(missing, ask = FALSE, update = FALSE)
  }

  # Load, and warn on version mismatches
  for (i in seq_along(pkgs)) {
    if (!requireNamespace(pkgs[i], quietly = TRUE)) {
      warning("Could not load '", pkgs[i], "'", call. = FALSE)
      next
    }
    if (check_version && !is.na(vers[i])) {
      installed <- as.character(packageVersion(pkgs[i]))
      if (installed != vers[i]) {
        warning(pkgs[i], ": installed ", installed, " but requirements.txt has ",
                vers[i], call. = FALSE)
      }
    }
    suppressPackageStartupMessages(
      library(pkgs[i], character.only = TRUE)
    )
  }
  invisible(pkgs)
}
