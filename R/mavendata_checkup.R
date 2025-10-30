#' mavendata_checkup
#'
#' @param file_path string path of mavendata file
#' @importFrom readxl read_excel
#' @importFrom tools file_path_sans_ext
#' @export
#'
mavendata_checkup <- function(file_path) {

  # Load file, skipping blank first row
  df_orig <- read_excel(file_path, skip = 1)
  messages <- list()
  setwd(dirname(file_path))
  # -------------------------------
  # QC PART 1: Check sample order
  # -------------------------------
  df <- df_orig[, -1]
  colnames_vector <- colnames(df)
  prefixes <- sub("_.*", "", colnames_vector)
  num_samples <- length(prefixes)
  lst <- list()

  for (j in 1:nrow(df)) {
    row_values <- as.character(unlist(df[j, ]))
    non_numeric <- sum(!grepl("^[0-9.]+$", row_values))
    has_sample_like <- any(grepl("QC|QS|neg|pos", row_values, ignore.case = TRUE))

    if (non_numeric > num_samples / 2 && has_sample_like) {
      row_prefixes <- sub("_.*", "", row_values)
      mismatch_idx <- which(row_prefixes != prefixes)

      if (length(mismatch_idx) > 0) {
        lst[[paste0("Row_", j)]] <- list(
          mismatches = mismatch_idx,
          expected = prefixes[mismatch_idx],
          found = row_prefixes[mismatch_idx]
        )
      }
    }
  }

  if (length(lst) == 0) {
    messages <- c(messages, "✅ No mismatches found in detected sample ID rows.")
  } else {
    messages <- c(messages, "❌ Mismatches found:")

    grouped <- list()
    make_key <- function(entry) {
      paste(
        paste(entry$mismatches, collapse = ","),
        paste(entry$expected, collapse = ","),
        paste(entry$found, collapse = ","),
        sep = "||"
      )
    }

    for (row_name in names(lst)) {
      entry <- lst[[row_name]]
      key <- make_key(entry)
      grouped[[key]] <- c(grouped[[key]], gsub("^Row_", "", row_name))
    }

    mismatch_messages <- character()
    for (key in names(grouped)) {
      parts <- strsplit(key, "\\|\\|")[[1]]
      mismatch_str <- parts[1]
      expected_str <- parts[2]
      found_str <- parts[3]
      rows_str <- paste(grouped[[key]], collapse = ", ")

      msg <- sprintf(
        "❌ Rows %s: Mismatches at positions [%s]\n  Expected: [%s]\n  Found:    [%s]",
        rows_str, mismatch_str, expected_str, found_str
      )
      mismatch_messages <- c(mismatch_messages, msg)
    }

    messages <- c(messages, paste(mismatch_messages, collapse = "\n\n"))
  }

  # -------------------------------
  # QC PART 2: Check spacing
  # -------------------------------
  check_column_pattern <- function(vec) {
    vec <- as.character(vec)
    vec[is.na(vec)] <- ""
    n <- length(vec)
    i <- 1
    block_num <- 1

    if (vec[i] != "") {
      messages <<- c(messages, "❌ Pattern violation: First row is not empty.")
      return(FALSE)
    }
    i <- i + 1

    while (i <= n) {
      if (i > n || vec[i] == "") {
        messages <<- c(messages, sprintf("❌ Pattern violation: Expected text at position %d.", i))
        return(FALSE)
      }
      while (i <= n && vec[i] != "") i <- i + 1

      empty_run <- 0
      while (i <= n && vec[i] == "") {
        empty_run <- empty_run + 1
        i <- i + 1
      }
      if (empty_run > 6) {
        messages <<- c(messages, sprintf("❌ End of data reached at row %d (more than 6 empty rows).", i - empty_run))
        return(TRUE)
      } else if (empty_run != 3 && i <= n) {
        messages <<- c(messages, sprintf("❌ Pattern violation: Expected 3 empty rows after text block %d, got %d.", block_num, empty_run))
        return(FALSE)
      }
      block_num <- block_num + 1
    }
    messages <<- c(messages, "✅ Column matches the expected pattern.")
    return(TRUE)
  }

  check_column_pattern(df_orig$...1)

  # -------------------------------
  # QC PART 3: Duplicate Metabolites
  # -------------------------------
  df_metab <- df_orig
  metabs <- na.omit(df_metab[, 1])
  names(metabs)[1] <- "Compound"
  dups <- unique(metabs$Compound[which(duplicated(metabs$Compound))])

  if (length(dups) > 0) {
    messages <- c(messages, sprintf("❌ There are duplicate metabolites of %s. Please fix.", dups))
    for (i in seq_along(dups)) {
      dup_rows <- na.omit(df_metab[which(grepl(dups[i], df_metab$...1)), ]) #FIX!!!!!!!
      if (any(duplicated(dup_rows))) {
        messages <- c(messages, sprintf("❌ Duplicates of %s have the same values", dups[i]))
      } else {
        messages <- c(messages, "❌ Duplicate metabs are not identical. Check peaks to determine correct values.")
      }
    }
  } else {
    messages <- c(messages, "✅ No duplicate metabolites")
  }

  # -------------------------------
  # QC PART 4: Extra RAW Files
  # -------------------------------

  raw_dir <- "RAW files"
  cat(file_path)
  # mzxml_dir <- "mzXML"
  extra_dir <- file.path(dirname(raw_dir), "extra_raw_files")

  if (!dir.exists(raw_dir)) stop("RAW files folder does not exist.")
  # if (!dir.exists(mzxml_dir)) stop("mzXML folder does not exist.")

  raw_files <- list.files(raw_dir, pattern = "\\.raw$", ignore.case = TRUE, full.names = FALSE)
  raw_basenames <- tools::file_path_sans_ext(raw_files)

  # ---- Code for if we want to run setdiff between raw files and unique mzxml files

  # mzxml_files <- list.files(mzxml_dir, pattern = "\\.mzXML$", ignore.case = TRUE, full.names = FALSE)
  # mzxml_basenames <- gsub("(_neg|_pos)?\\.mzXML$", "", mzxml_files, ignore.case = TRUE)
  # mzxml_unique <- unique(mzxml_basenames)
  #
  # extra_raw_basenames <- setdiff(raw_basenames, mzxml_unique)

  #extra_raw_basenames <- setdiff(raw_basenames, prefixes)
  # Identify overlap
  common <- intersect(raw_basenames, prefixes)

  if (length(common) == 0) {
    messages <- c(messages, sprintf("⚠️ No overlap found between raw_basenames and prefixes."))
    extra_raw_basenames <- character(0)  # or raw_basenames if you want to treat all as extras
  } else {
    extra_raw_basenames <- raw_basenames[!(raw_basenames %in% prefixes)]
  }
  extra_raw_fullpaths <- file.path(raw_dir, paste0(extra_raw_basenames, ".raw"))

  if (length(extra_raw_basenames) > 0) {
    if (!dir.exists(extra_dir)) dir.create(extra_dir)
    for (f in extra_raw_fullpaths) {
      file.rename(from = f, to = file.path(extra_dir, basename(f)))
    }
    messages <- c(
      messages,
      paste0(
        "❌ but ✅ Moved the following extra RAW files to extra_raw_files:\n  ",
        paste(basename(extra_raw_fullpaths), collapse = "\n  ")
      )
    )
  } else {
    messages <- c(messages, "✅ No extra RAW files found.")
  }

  # -------------------------------
  # Final Report
  # -------------------------------
  cat("\n=========== Mavendata Checkup Report ===========\n\n")
  for (msg in messages) {
    cat(paste(msg, collapse = "\n"), "\n\n")
  }
  cat("======================================\n")
}
