# Manuscript

The manuscript is `main.tex`, with references in `references.bib` and author details in `author.tex`.
The six current figures use PNG assets in `figures/`.
Conceptual figure sources are retained as LaTeX files.

Build from this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

See the repository README for commands to regenerate the result figures.
Check author information, repository accessibility, data permissions, and journal requirements before submission.
Local builds do not synchronize Overleaf.
