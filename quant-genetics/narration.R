library(stringr)
library(fs)
library(jsonlite)

# =========================
# CONFIGURATION
# =========================

slides_dir <- "slides"
narr_dir   <- "narration"

voice_model <- "gpt-4o-mini-tts"
voice_name  <- "alloy"
update_all  <- FALSE   # TRUE = regenerate audio

dir_create(narr_dir, recurse = TRUE)

api_key <- Sys.getenv("OPENAI_API_KEY")
if (api_key == "") stop("OPENAI_API_KEY not set.")

qmd_files <- dir_ls(slides_dir, regexp = "\\.qmd$")
qmd_files <- qmd_files[!str_detect(qmd_files, "_narrated\\.qmd$")]

# Only first file (testing mode)
qmd_files <- qmd_files[1]

# =========================
# HELPERS
# =========================

clean_text_for_tts <- function(text) {
  text |>
    str_replace_all("\\$[^$]+\\$", "") |>
    str_replace_all("\\\\[a-zA-Z]+", "") |>
    str_replace_all("\\{.*?\\}", "") |>
    str_squish()
}

generate_audio <- function(text, output_path) {
  payload <- toJSON(
    list(
      model = voice_model,
      voice = voice_name,
      input = text
    ),
    auto_unbox = TRUE
  )
  
  cmd <- paste(
    "curl https://api.openai.com/v1/audio/speech",
    "-H", shQuote(paste0("Authorization: Bearer ", api_key)),
    "-H", shQuote("Content-Type: application/json"),
    "-d", shQuote(payload),
    "--output", shQuote(output_path)
  )
  
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  file_exists(output_path) && file_info(output_path)$size > 1000
}

remove_existing_audio_tags <- function(lines) {
  lines[!str_detect(lines, "narration-player|<audio\\b")]
}

insert_audio_tag_before_notes <- function(lines, note_start, filename) {
  
  audio_tag <- sprintf(
    '<audio controls preload="auto" src="../%s/%s"></audio>',
    narr_dir, filename
  )
  
  append(lines, c("", audio_tag, ""), after = note_start - 1)
}

render_slide_copy <- function(narrated_qmd, out_dir = slides_dir) {
  narrated_qmd <- normalizePath(narrated_qmd, winslash = "/", mustWork = TRUE)
  dest_dir     <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  
  tmp_proj <- "slides_render"
  dir.create(tmp_proj, showWarnings = FALSE)
  
  yml <- file.path(tmp_proj, "_quarto.yml")
  if (!file.exists(yml)) {
    writeLines(c("project:", "  type: default"), yml)
  }
  
  tmp_qmd <- file.path(tmp_proj, basename(narrated_qmd))
  file.copy(narrated_qmd, tmp_qmd, overwrite = TRUE)
  
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(tmp_proj)
  
  res <- system2("quarto", c("render", basename(tmp_qmd)))
  if (res != 0) stop("Render failed for: ", narrated_qmd)
  
  html <- sub("\\.qmd$", ".html", basename(tmp_qmd))
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  file.copy(html, file.path(dest_dir, html), overwrite = TRUE)
  
  invisible(file.path(dest_dir, html))
}

# =========================
# MAIN PIPELINE
# =========================

summary_report <- list()

for (qfile in qmd_files) {
  
  message("\nProcessing: ", qfile)
  
  original_lines <- readLines(qfile, warn = FALSE)
  original_lines <- remove_existing_audio_tags(original_lines)
  
  new_lines <- original_lines
  
  note_starts <- which(str_detect(new_lines, "^::: notes"))
  base_name   <- path_ext_remove(path_file(qfile))
  
  if (length(note_starts) == 0) {
    message("  No notes found.")
    next
  }
  
  blocks_generated <- 0
  
  # 🔥 PROCESS FROM BOTTOM TO TOP (so line numbers stay valid)
  for (i in rev(seq_along(note_starts))) {
    
    start <- note_starts[i]
    
    # Find matching closing ::: for text extraction only
    end <- NA
    for (j in (start + 1):length(new_lines)) {
      if (str_detect(new_lines[j], "^:::\\s*$")) {
        end <- j
        break
      }
    }
    if (is.na(end)) next
    
    raw_text  <- paste(new_lines[(start + 1):(end - 1)], collapse = " ")
    clean_txt <- clean_text_for_tts(raw_text)
    
    out_name <- sprintf("%s_block%02d.mp3", base_name, i)
    out_path <- path(narr_dir, out_name)
    
    if (!file_exists(out_path) || update_all) {
      message("  Generating audio block ", i)
      success <- generate_audio(clean_txt, out_path)
      if (!success) {
        warning("  Failed block ", i, " (", out_name, ")")
        next
      }
      blocks_generated <- blocks_generated + 1
    } else {
      message("  Using existing audio block ", i)
    }
    
    # ✅ INSERT BEFORE notes block
    new_lines <- insert_audio_tag_before_notes(
      new_lines,
      note_start = start,
      filename = out_name
    )
  }
  
  narrated_qmd <- path(slides_dir, paste0(base_name, "_narrated.qmd"))
  writeLines(new_lines, narrated_qmd)
  
  message("  Rendering narrated HTML...")
  render_slide_copy(narrated_qmd, out_dir = slides_dir)
  
  summary_report[[base_name]] <- list(
    blocks    = length(note_starts),
    generated = blocks_generated,
    output    = paste0(base_name, "_narrated.html")
  )
}

# =========================
# SUMMARY
# =========================

cat("\n====================================\n")
cat("  NARRATION PIPELINE SUMMARY\n")
cat("====================================\n\n")

for (name in names(summary_report)) {
  cat("Slide deck:", name, "\n")
  cat("  Notes blocks:    ", summary_report[[name]]$blocks, "\n")
  cat("  Audio generated: ", summary_report[[name]]$generated, "\n")
  cat("  Output HTML:     ", file.path(slides_dir, summary_report[[name]]$output), "\n\n")
}

cat("Pipeline complete.\n")

# library(stringr)
# library(fs)
# library(jsonlite)
# library(quarto)
# 
# #setwd("C://Users//au223366//Documents//GitHub//qgteach//quant-genetics")
# #getwd()
# 
# # =========================
# # CONFIGURATION
# # =========================
# 
# slides_dir <- "slides"
# narr_dir   <- "narration"
# 
# voice_model <- "gpt-4o-mini-tts"
# voice_name  <- "alloy"
# update_all  <- FALSE   # TRUE = regenerate audio
# 
# dir_create(narr_dir, recurse = TRUE)
# 
# api_key <- Sys.getenv("OPENAI_API_KEY")
# if (api_key == "") stop("OPENAI_API_KEY not set.")
# 
# qmd_files <- dir_ls(slides_dir, regexp = "\\.qmd$")
# qmd_files <- qmd_files[!str_detect(qmd_files, "_narrated\\.qmd$")]
# 
# qmd_files <- qmd_files[1]
# 
# # =========================
# # HELPERS
# # =========================
# 
# clean_text_for_tts <- function(text) {
#   text |>
#     str_replace_all("\\$[^$]+\\$", "") |>     # remove inline math
#     str_replace_all("\\\\[a-zA-Z]+", "") |>   # remove LaTeX commands
#     str_replace_all("\\{.*?\\}", "") |>       # remove braces
#     str_squish()
# }
# 
# generate_audio <- function(text, output_path) {
#   
#   payload <- toJSON(
#     list(
#       model = voice_model,
#       voice = voice_name,
#       input = text
#     ),
#     auto_unbox = TRUE
#   )
#   
#   cmd <- paste(
#     "curl https://api.openai.com/v1/audio/speech",
#     "-H", shQuote(paste0("Authorization: Bearer ", api_key)),
#     "-H", shQuote("Content-Type: application/json"),
#     "-d", shQuote(payload),
#     "--output", shQuote(output_path)
#   )
#   
#   system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
#   
#   file_exists(output_path) && file_info(output_path)$size > 1000
# }
# 
# #remove_existing_audio_tags <- function(lines) {
# #  lines[!str_detect(lines, "<audio controls")]
# #}
# 
# remove_existing_audio_tags <- function(lines) {
#   lines[!str_detect(lines, "narration-player|<audio\\b")]
# }
# 
# # insert_audio_tag <- function(lines, insert_after, filename) {
# #   
# #   audio_tag <- sprintf(
# #     '<audio controls style="width:280px;"><source src="../%s/%s" type="audio/mpeg"></audio>',
# #     narr_dir,
# #     filename
# #   )
# #   
# #   append(lines, audio_tag, after = insert_after)
# # }
# 
# # insert_audio_tag <- function(lines, insert_after, filename) {
# #   audio_tag <- sprintf(
# #     '<div class="narration-player"><audio controls><source src="../%s/%s" type="audio/mpeg"></audio></div>',
# #     narr_dir,
# #     filename
# #   )
# #   append(lines, audio_tag, after = insert_after)
# # }
# 
# insert_audio_tag <- function(lines, note_start, note_end, filename) {
#   audio_tag <- sprintf(
#     '<div class="narration-player"><audio controls><source src="../%s/%s" type="audio/mpeg"></audio></div>',
#     narr_dir, filename
#   )
#   
#   # Insert AFTER the closing ::: of the notes block (safe place)
#   block <- c("", audio_tag, "")
#   append(lines, block, after = note_end)
# }
# 
# render_slide_copy <- function(narrated_qmd, out_dir = "slides") {
#   narrated_qmd <- normalizePath(narrated_qmd, winslash = "/", mustWork = TRUE)
#   
#   tmp_proj <- "slides_render"
#   dir.create(tmp_proj, showWarnings = FALSE)
#   
#   # Ensure neutral _quarto.yml exists
#   yml <- file.path(tmp_proj, "_quarto.yml")
#   if (!file.exists(yml)) {
#     writeLines(c("project:", "  type: default"), yml)
#   }
#   
#   # Copy qmd into neutral project
#   tmp_qmd <- file.path(tmp_proj, basename(narrated_qmd))
#   file.copy(narrated_qmd, tmp_qmd, overwrite = TRUE)
#   
#   # Render in neutral project
#   old <- getwd()
#   on.exit(setwd(old), add = TRUE)
#   setwd(tmp_proj)
#   
#   res <- system2("quarto", c("render", basename(tmp_qmd)))
#   if (res != 0) stop("Render failed")
#   
#   # Move HTML back to slides/
#   html <- sub("\\.qmd$", ".html", basename(tmp_qmd))
#   file.copy(html, file.path(out_dir, html), overwrite = TRUE)
#   
#   invisible(file.path(out_dir, html))
# }
# 
# 
# # =========================
# # MAIN PIPELINE
# # =========================
# 
# summary_report <- list()
# 
# for (qfile in qmd_files) {
#   
#   message("\nProcessing: ", qfile)
#   
#   original_lines <- readLines(qfile, warn = FALSE)
#   original_lines <- remove_existing_audio_tags(original_lines)
#   
#   new_lines <- original_lines
#   
#   note_starts <- which(str_detect(new_lines, "^::: notes"))
#   base_name   <- path_ext_remove(path_file(qfile))
#   
#   if (length(note_starts) == 0) {
#     message("  No notes found.")
#     next
#   }
#   
#   blocks_generated <- 0
#   offset <- 0
#   
#   for (i in seq_along(note_starts)) {
#     
#     start <- note_starts[i] + offset
#     
#     # Find closing :::
#     end <- NA
#     for (j in (start + 1):length(new_lines)) {
#       if (str_detect(new_lines[j], "^:::\\s*$")) {
#         end <- j
#         break
#       }
#     }
#     
#     if (is.na(end)) next
#     
#     raw_text  <- paste(new_lines[(start + 1):(end - 1)], collapse = " ")
#     clean_txt <- clean_text_for_tts(raw_text)
#     
#     out_name <- sprintf("%s_block%02d.mp3", base_name, i)
#     out_path <- path(narr_dir, out_name)
#     
#     if (!file_exists(out_path) || update_all) {
#       message("  Generating audio block ", i)
#       success <- generate_audio(clean_txt, out_path)
#       if (!success) {
#         warning("  Failed block ", i)
#         next
#       }
#       blocks_generated <- blocks_generated + 1
#     } else {
#       message("  Using existing audio block ", i)
#     }
#     
#     # Insert audio tag before ::: notes
#     #new_lines <- insert_audio_tag(new_lines, start - 1, out_name)
#     #offset <- offset + 1
#     # Insert audio tag AFTER the ::: notes block (safe place)
#     new_lines <- insert_audio_tag(new_lines, note_start = start, note_end = end, filename = out_name)
#     
#     # We inserted 3 lines: "", audio_tag, ""
#     offset <- offset + 3
#   }
#   
#   narrated_qmd <- path(slides_dir, paste0(base_name, "_narrated.qmd"))
#   
#   writeLines(new_lines, narrated_qmd)
#   
#   message("  Rendering narrated HTML...")
#   #quarto_render(narrated_qmd)
#   #quarto_render(
#   #  input = narrated_qmd,
#   #  project = FALSE
#   #)
#   
#   #res <- system2("quarto", c("render", narrated_qmd, "--no-project"))
#   #if (res != 0) warning("Quarto render failed for: ", narrated_qmd)
#   
# 
#   #render_slide_copy("slides/introduction_population_genetics_narrated.qmd")
#   render_slide_copy(narrated_qmd, out_dir = slides_dir)
# 
#   summary_report[[base_name]] <- list(
#     blocks = length(note_starts),
#     generated = blocks_generated,
#     output = paste0(base_name, "_narrated.html")
#   )
# }
# 
# # =========================
# # SUMMARY
# # =========================
# 
# cat("\n====================================\n")
# cat("  NARRATION PIPELINE SUMMARY\n")
# cat("====================================\n\n")
# 
# for (name in names(summary_report)) {
#   cat("Slide deck:", name, "\n")
#   cat("  Notes blocks:    ", summary_report[[name]]$blocks, "\n")
#   cat("  Audio generated: ", summary_report[[name]]$generated, "\n")
#   cat("  Output HTML:     slides/", summary_report[[name]]$output, "\n\n")
# }
# 
# cat("Pipeline complete.\n")