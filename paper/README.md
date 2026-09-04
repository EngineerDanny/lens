# Manuscript source

The archive-ready LaTeX source is `main.tex`. The bibliography is in `references.bib`, and author metadata is isolated in `author.tex`.

Compile with:

```sh
latexmk -pdf main.tex
```

Before public submission, replace `Anonymous manuscript` in `author.tex` with the final author and affiliation block, and add the permanent public code repository address in the Data and code availability section. No scientific text, results, figures, or references depend on those two metadata changes.

