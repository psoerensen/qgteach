# qgteach

Teaching materials for quantitative genetics, genomics, and statistical modeling using R.

This repository contains course materials developed using **R** and **Quarto**, including:

- lecture slides
- theoretical notes
- computational tutorials
- interactive demonstrations
- datasets and figures

The materials emphasize **reproducible research**, **statistical modeling**, and **computational implementation** in quantitative genetics and genomics.

Website  
https://psoerensen.github.io/qgteach/


## Repository Structure

The repository is organized as a teaching website containing multiple courses and shared tools.

```
qgteach/

  _quarto.yml
  index.qmd
  about.qmd

  docs/                      # rendered website (GitHub Pages)

  tools/                     # utilities for course creation and narration
    create_course.R
    narrate_slides_qmd.R
    styles.css

  quant-genetics/            # course materials
  genomics-systems-bioinfo/  # course materials
```

Each course has its own folder and Quarto configuration.


## Course Structure

Each course typically follows this structure:

```
course-name/

  _quarto.yml
  index.qmd
  styles.css

  slides/        # lecture slides
  notes/         # theoretical notes
  tutorials/     # hands-on exercises
  apps/          # interactive teaching apps

  data/          # datasets used in tutorials
  images/        # figures used in slides and notes

  narration/     # narration audio files
```

Slides can optionally include speaker notes that are used for generating narrated lectures.


## Teaching Philosophy

The teaching material emphasizes:

- mathematical foundations of quantitative genetics
- statistical modeling
- computational implementation in R
- reproducible research workflows
- conceptual understanding of genomic analyses

Lecture materials integrate **theory**, **simulation**, and **applied data analysis**.


## Software Used

The courses make use of several R packages developed for genomic analysis.

### qgg

https://github.com/psoerensen/qgg

Provides tools for:

- Bayesian linear regression models
- genomic fine mapping
- polygenic scoring
- gene set enrichment analysis


### gact

https://github.com/psoerensen/gact

Provides tools for:

- genomic association databases
- GWAS summary statistics processing
- linking genetic variants to biological pathways
- integrative genomics workflows


## Rendering the Website

The teaching site is built using **Quarto**.

To render the full website locally:

```
quarto render
```

The rendered site is written to:

```
docs/
```

This folder is used for **GitHub Pages hosting** and should normally not be edited manually.


## Creating a New Course

New courses can be created using the helper function located in `tools/create_course.R`.

From the repository root in R:

```r
source("tools/create_course.R")
create_course("my-new-course")
```

This automatically creates the course structure:

```
my-new-course/

  _quarto.yml
  index.qmd
  styles.css

  slides/
  notes/
  tutorials/
  apps/
  data/
  images/
  narration/
```

It also creates example files:

```
slides/01_introduction.qmd
tutorials/01_getting_started.qmd
```

The slide template includes a **notes block** that can be used for narration.


## Generating Narrated Slides

Slides can include narration using Reveal.js notes blocks:

```
:::
Explain the motivation for this topic.
Provide intuition and background.
Mention key points students should remember.
:::
```

Narrated audio can be generated automatically using the script in:

```
tools/narrate_slides_qmd.R
```

Example usage:

```r
source("tools/narrate_slides_qmd.R")

narrate_slides("quant-genetics")
```

The script will:

- extract text from slide notes
- generate audio using the OpenAI TTS API
- save audio files in the `narration/` folder
- create narrated slide versions


## Contributing

Contributions are welcome.

Possible contributions include:

- new courses
- new lectures
- new tutorials
- improved explanations
- additional datasets
- bug fixes

When adding new material:

- slides → `slides/`
- theory notes → `notes/`
- exercises → `tutorials/`
- small datasets → `data/`

Please keep tutorials **reproducible** and avoid committing large datasets.


## Maintainer

Peter Sørensen  
Center for Quantitative Genetics and Genomics  
Aarhus University  

Email: pso@qgg.au.dk  
GitHub  
https://github.com/psoerensen
