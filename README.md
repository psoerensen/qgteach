# qgteach

![Quarto](https://img.shields.io/badge/built%20with-Quarto-blue)
![R](https://img.shields.io/badge/language-R-blue)

Teaching materials for **quantitative genetics, genomics, and statistical modeling using R**.

This repository contains modular course materials developed using **R** and **Quarto**, including:

- lecture slides
- theoretical notes
- tutorials and exercises
- interactive demonstrations

The repository also provides tools for creating new courses and generating narrated lecture slides.

Website  
https://psoerensen.github.io/qgteach/

## Repository Structure

```text
qgteach/

  _quarto.yml
  index.qmd
  about.qmd

  docs/                      # rendered website (GitHub Pages)

  tools/                     # course scaffolding and narration utilities
    create_course.R
    narrate_slides_qmd.R
    styles.css

  quant-genetics/
  genomics-systems-bioinfo/
```

## Course Structure

Each course follows a standard layout.

```text
course-name/

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

Slides may include **speaker notes** used to generate narrated lectures.

## Tools

### create_course.R

Creates a new course with the standard structure and template files.

```r
source("tools/create_course.R")
create_course("my-course")
```

Templates include:

```text
slides/01_introduction.qmd
notes/01_concepts.qmd
tutorials/01_getting_started.qmd
exercises/01_exercise_template.qmd
```

### narrate_slides_qmd.R

Generates narrated slides from Reveal.js notes blocks.

Example notes block:

```text
:::
Explain the key idea for this slide.
:::
```

Generate narration:

```r
source("tools/narrate_slides_qmd.R")
narrate_slides("quant-genetics")
```

Audio files are saved to the course `narration/` folder.

## Rendering the Website

Render the full teaching site with:

```bash
quarto render
```

The output is written to:

```text
docs/
```

which is used for **GitHub Pages hosting**.

## Typical Workflow

1. Create a new course

```r
create_course("my-course")
```

2. Write slides, notes, tutorials, and exercises

3. Add narration notes to slides

4. Generate narrated slides

```r
narrate_slides("my-course")
```

5. Render the site

```bash
quarto render
```

## Teaching Approach

The material emphasizes:

- mathematical foundations of quantitative genetics
- statistical modeling
- reproducible workflows in R
- integration of theory, simulation, and real data

## Contributing

Contributions are welcome.

Examples include:

- new courses
- lectures or tutorials
- exercises
- improved explanations
- bug fixes

Please keep tutorials and exercises **reproducible** and avoid committing large datasets.

## Maintainer

Peter Sørensen  
Center for Quantitative Genetics and Genomics  
Aarhus University  

Email: pso@qgg.au.dk  
GitHub  
https://github.com/psoerensen
