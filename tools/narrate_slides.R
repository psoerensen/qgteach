narrate_slides <- function(course_dir = ".",
                           voice_model = "gpt-4o-mini-tts",
                           voice_name  = "alloy",
                           update_all  = FALSE) {
  
  library(stringr)
  library(fs)
  library(jsonlite)
  
  slides_dir <- file.path(course_dir, "slides")
  narr_dir   <- file.path(course_dir, "narration")
  
  dir_create(narr_dir, recurse = TRUE)
  
  api_key <- Sys.getenv("OPENAI_API_KEY")
  if (api_key == "") stop("OPENAI_API_KEY not set.")
  
  qmd_files <- dir_ls(slides_dir, regexp = "\\.qmd$")
  qmd_files <- qmd_files[!str_detect(qmd_files, "_narrated\\.qmd$")]
  
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
      '<audio controls preload="auto" src="../narration/%s"></audio>',
      filename
    )
    
    append(lines, c("", audio_tag, ""), after = note_start - 1)
  }
  
  # =========================
  # MAIN LOOP
  # =========================
  
  for (qmd in qmd_files) {
    
    message("Processing: ", basename(qmd))
    
    lines <- readLines(qmd)
    
    # remove old audio
    lines <- remove_existing_audio_tags(lines)
    
    text <- paste(lines, collapse = " ")
    text <- clean_text_for_tts(text)
    
    audio_file <- file.path(
      narr_dir,
      paste0(tools::file_path_sans_ext(basename(qmd)), ".mp3")
    )
    
    if (!file_exists(audio_file) || update_all) {
      
      ok <- generate_audio(text, audio_file)
      
      if (!ok) {
        warning("Audio generation failed: ", qmd)
        next
      }
      
    }
    
    note_line <- which(str_detect(lines, ":::\\s*notes"))[1]
    
    if (!is.na(note_line)) {
      
      filename <- basename(audio_file)
      
      lines <- insert_audio_tag_before_notes(lines, note_line, filename)
      
      narrated_qmd <- sub("\\.qmd$", "_narrated.qmd", qmd)
      
      writeLines(lines, narrated_qmd)
      
      message("Created narrated slide: ", narrated_qmd)
      
    }
    
  }
  
  message("Narration generation complete.")
}