# qgteach

![Quarto](https://img.shields.io/badge/built%20with-Quarto-blue)
![R](https://img.shields.io/badge/language-R-blue)

Teaching materials for quantitative genetics, genomics, and statistical modeling using R.

This repository contains modular course materials developed using **R** and **Quarto**, including lecture slides, theoretical notes, computational tutorials, exercises, and interactive demonstrations.

The repository also provides tools for creating new courses and generating narrated lecture slides.

Website  
https://psoerensen.github.io/qgteach/


## Repository Structure

The repository is organized as a Quarto website containing multiple courses and shared tools.

```
qgteach/

  _quarto.yml
  index.qmd
  about.qmd

  docs/                      # rendered website (GitHub Pages)

  tools/                     # utilities for course creation and narration
    create_course.R
    narrate_slides.R
    styles.css

  quant-genetics/            # course materials
  genomics-systems-bioinfo/  # course materials
```

Each course has its own folder and Quarto configuration.


## Course Structure

Each course follows a standard structure.

```
course-name/

  _quarto.yml
  index.qmd
  styles.css

  slides/        # lecture slides
  notes/         # theoretical notes

  tutorials/     # guided practical walkthroughs
  exercises/     # student problem sets

  apps/          # interactive teaching apps

  data/          # datasets used in tutorials and exercises
  images/        # figures used in slides and notes

  narration/     # narration audio files
```

Slides may optionally contain **speaker notes** that can be used to generate narrated lectures.


## Creating a New Course

New courses can be created using the helper function in `tools/create_course.R`.

Example:

```r
source("tools/create_course.R")
create_course("my-new-course")
```

This automatically generates the course structure:

```
my-new-course/

  _quarto.yml
  index.qmd
  styles.css

  slides/
  notes/
  tutorials/
  exercises/
  apps/
  data/
  images/
  narration/
```

Template files are also created:

```
slides/01_introduction.qmd
tutorials/01_getting_started.qmd
exercises/01_exercise_template.qmd
```

The slide template includes a **notes block** that can be used for narration.


## Generating Narrated Slides

Slides can include narration using Reveal.js notes blocks:

```
:::
Explain the motivation for this topic.
Provide intuition and background.
Mention key ideas students should remember.
:::
```

Narrated audio can then be generated using:

```r
source("tools/narrate_slides.R")
narrate_slides("quant-genetics")
```

The script will:

- extract text from slide notes
- generate audio using the OpenAI TTS API
- save audio files in the `narration/` folder
- create narrated slide versions


## Rendering the Website

The teaching site is built using **Quarto**.

To render the entire site locally:

```
quarto render
```

The rendered website is written to:

```
docs/
```

This directory is used for **GitHub Pages hosting**.


## Teaching Approach

The material emphasizes:

- mathematical foundations of quantitative genetics
- statistical modeling
- reproducible computational workflows in R
- integration of theory, simulation, and real data

Courses combine conceptual explanations with interactive examples and computational exercises.


## Contributing

Contributions are welcome.

Examples of contributions include:

- new courses
- additional lectures
- new tutorials
- new exercises
- improved explanations
- bug fixes

When adding new material:

- slides → `slides/`
- theory notes → `notes/`
- guided walkthroughs → `tutorials/`
- problem sets → `exercises/`
- small datasets → `data/`

Please keep tutorials and exercises **reproducible** and avoid committing large datasets.


## Maintainer

Peter Sørensen  
Center for Quantitative Genetics and Genomics  
Aarhus University  

Email: pso@qgg.au.dk  
GitHub  
https://github.com/psoerensen
