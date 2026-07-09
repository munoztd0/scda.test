# 2. Identify Target Scripts
# 'ls*.R' matches any file starting with 'ls' and ending in '.R'
all_scripts <- list.files(path = "source_code/", pattern = "^lsfvit.*\\.R$", full.names = TRUE)
#tsvit #tsi
#  3. Execution Loop
for (script in all_scripts) {
  message("Processing: ", script)

  # SAFETY: Clear the specific object from the previous run
  if (exists("result")) rm(result)
  result <- NULL

  # 4. Attempt to process the file
  tryCatch({
    # Run the script in a local environment to find 'result'
    # we use a new environment to prevent global variable pollution
    script_env <- new.env()
    source(script, local = script_env)

    # Extract 'result' from that script's environment
    if (!is.null(script_env$result)) {

      # Calculate widths
      new_widths <- junco::def_colwidths(head(script_env$result,100), fontspec = font_spec("Times", 9L, 1.2), label_width_ins = 2)



      # Format into R code string
      code_string <- paste0("colwidth <- ", paste(deparse(new_widths), collapse = ""))

      # Read and Replace
      lines <- readLines(script)
      target_idx <- grep("# \\[AUTO-COLWIDTH\\]", lines)

      if (length(target_idx) > 0) {
        lines[target_idx] <- code_string
        writeLines(lines, script)
        message("  -> Success: Injected ", code_string)
      } else {
        message("  -> Warning: No placeholder found in ", script)
      }

    } else {
      message("  -> Error: 'result' object not found after sourcing ", script)
    }

  }, error = function(e) {
    stop("  -> Critical Failure in ", script, ": ", e$message)
  })

  # Cleanup for next iteration
  rm(script_env)
}

