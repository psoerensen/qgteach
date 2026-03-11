lines <- readLines("slides/introduction_population_genetics.qmd", warn = FALSE)

open_stack <- character()
problems <- list()

is_open <- function(x) grepl("^:::\\s*\\{", x) || grepl("^::: notes\\s*$", x)
is_close <- function(x) grepl("^:::\\s*$", x)  # EXACT close only

for (i in seq_along(lines)) {
  x <- lines[i]
  
  if (is_open(x)) {
    open_stack <- c(open_stack, paste0(i, ": ", x))
  } else if (is_close(x)) {
    if (length(open_stack) == 0) {
      problems <- c(problems, list(paste("Stray close at line", i)))
    } else {
      open_stack <- open_stack[-length(open_stack)]
    }
  } else if (grepl("^:::\\s+$", x)) {
    problems <- c(problems, list(paste("Close fence has trailing spaces at line", i, "->", shQuote(x))))
  }
}

cat("---- Problems ----\n")
if (length(problems)) cat(paste0("* ", unlist(problems), collapse="\n"), "\n") else cat("None found.\n")

cat("\n---- Unclosed openings left on stack ----\n")
if (length(open_stack)) cat(paste0("* ", open_stack, collapse="\n"), "\n") else cat("None.\n")
