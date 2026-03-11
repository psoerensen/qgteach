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

The repository is organized as a teaching website containing multiple courses.

```
qgteach/

  _quarto.yml
  index.qmd
  about.qmd

  docs/                      # rendered website (GitHub Pages)

  quant-genetics/            # course materials
  genomics-systems-bioinfo/  # course materials

  narration/                 # shared narration tools
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

  narration/     # narration files for lectures
```


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

To render locally:

```
quarto render
```

The rendered site is written to:

```
docs/
```

This folder is used for **GitHub Pages hosting** and should normally not be edited manually.


## Contributing

Contributions are welcome.

Possible contributions include:

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


## Adding a New Course

To add a new course to the repository:

1. Create a new folder in the repository root.
2. Add a course-specific `_quarto.yml`.
3. Add an `index.qmd` file for the course homepage.
4. Follow the standard course structure.

Example:

```
new-course/

  _quarto.yml
  index.qmd

  slides/
  notes/
  tutorials/
  apps/
  data/
  images/
```

After adding the course, update the main `index.qmd` page to link to it.


## Maintainer

Peter Sørensen  
Center for Quantitative Genetics and Genomics  
Aarhus University  

Email: pso@qgg.au.dk  
GitHub  
https://github.com/psoerensen