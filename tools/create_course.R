create_course <- function(name, path = ".", tools_dir = "tools") {
  
  course_dir <- file.path(path, name)
  
  if (!dir.exists(course_dir)) {
    dir.create(course_dir, recursive = TRUE)
    message("Created course directory: ", course_dir)
  }
  
  dirs <- c(
    "slides",
    "notes",
    "tutorials",
    "exercises",
    "apps",
    "data",
    "images",
    "narration"
  )
  
  for (d in dirs) {
    
    dir_path <- file.path(course_dir, d)
    
    if (!dir.exists(dir_path)) {
      dir.create(dir_path)
      message("Created folder: ", dir_path)
    }
    
  }
  
  # ----------------------------
  # index.qmd
  # ----------------------------
  
  index_file <- file.path(course_dir, "index.qmd")
  
  if (!file.exists(index_file)) {
    
    index <- c(
      "---",
      paste0('title: "', name, '"'),
      "---",
      "",
      "## Course Material",
      "",
      "Add a short description of the course.",
      "",
      "## Slides",
      "",
      "Slides are located in the `slides/` folder.",
      "",
      "## Tutorials",
      "",
      "Guided tutorials are located in the `tutorials/` folder.",
      "",
      "## Exercises",
      "",
      "Student problem sets are located in the `exercises/` folder.",
      "",
      "## Interactive Tools",
      "",
      "Apps and simulations are located in the `apps/` folder.",
      "",
      "## Prerequisites",
      "",
      "- Basic genetics",
      "- Introductory statistics",
      "- Familiarity with R"
    )
    
    writeLines(index, index_file)
    
    message("Created file: ", index_file)
    
  }
  
  # ----------------------------
  # _quarto.yml
  # ----------------------------
  
  quarto_file <- file.path(course_dir, "_quarto.yml")
  
  if (!file.exists(quarto_file)) {
    
    quarto <- c(
      "project:",
      "  type: book",
      "",
      "resources:",
      "  - images/**",
      "  - narration/**",
      "  - slides/**",
      "  - apps/**",
      "  - data/**",
      "",
      "book:",
      paste0('  title: "', name, '"'),
      "  chapters:",
      "    - index.qmd",
      "",
      "format:",
      "  html:",
      "    theme: cosmo",
      "    toc: true"
    )
    
    writeLines(quarto, quarto_file)
    
    message("Created file: ", quarto_file)
    
  }
  
  # ----------------------------
  # copy styles.css
  # ----------------------------
  
  css_file <- file.path(course_dir, "styles.css")
  
  if (!file.exists(css_file)) {
    
    source_css <- file.path(tools_dir, "styles.css")
    
    if (file.exists(source_css)) {
      
      file.copy(source_css, css_file)
      
      message("Copied styles.css")
      
    } else {
      
      warning("styles.css not found in tools directory")
      
    }
    
  }
  
  # ----------------------------
  # template slide
  # ----------------------------
  
  slide_file <- file.path(course_dir, "slides", "01_introduction.qmd")
  
  if (!file.exists(slide_file)) {
    
    slide <- c(
      "---",
      'title: "Introduction"',
      "format: revealjs",
      "---",
      "",
      "# Introduction",
      "",
      "Add lecture content here.",
      "",
      "::: notes",
      "Explain the motivation for this topic.",
      "",
      "Provide background context and intuition.",
      "",
      "Mention key ideas students should focus on.",
      ":::"
    )
    
    writeLines(slide, slide_file)
    
    message("Created template slide")
    
  }
  
  # ----------------------------
  # template tutorial
  # ----------------------------
  
  tutorial_file <- file.path(course_dir, "tutorials", "01_getting_started.qmd")
  
  if (!file.exists(tutorial_file)) {
    
    tutorial <- c(
      "---",
      'title: "Getting Started"',
      "format: html",
      "---",
      "",
      "## Tutorial",
      "",
      "Example R code:",
      "",
      "```r",
      "x <- rnorm(100)",
      "mean(x)",
      "```"
    )
    
    writeLines(tutorial, tutorial_file)
    
    message("Created template tutorial")
    
  }
  
  # ----------------------------
  # template exercise
  # ----------------------------
  
  exercise_file <- file.path(course_dir, "exercises", "01_exercise_template.qmd")
  
  if (!file.exists(exercise_file)) {
    
    exercise <- c(
      "---",
      'title: "Exercise"',
      "format: html",
      "---",
      "",
      "## Exercise",
      "",
      "Describe the problem students should solve.",
      "",
      "### Tasks",
      "",
      "1. Implement the model in R.",
      "2. Visualize the results.",
      "3. Interpret the output.",
      "",
      "```r",
      "# write your code here",
      "```"
    )
    
    writeLines(exercise, exercise_file)
    
    message("Created template exercise")
    
  }
  
  message("Course setup complete.")
}