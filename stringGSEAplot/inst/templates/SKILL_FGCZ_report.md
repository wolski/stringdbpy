# FGCZ Quarto Report Style Guide

## Purpose

This document defines the FGCZ report conventions for Quarto `.qmd` files.
Follow these patterns when creating new FGCZ-branded reports.

Reference implementation: `ezRun/inst/templates/quarto/diffExp_proteomics.qmd`
and `diffExp_rnaseq.qmd`.

---

## Shared Infrastructure

Every FGCZ Quarto report uses these shared files (in the same directory as the `.qmd`):

| File | Purpose |
|---|---|
| `_fgcz-report.yml` | Shared Quarto metadata: self-contained HTML, cosmo theme, FGCZ header, lightbox, knitr defaults, code-fold |
| `fgcz_header_quarto.html` | FGCZ banner (base64 image), inline CSS for banner/tabs/callouts, Enrichr JS, external-link JS |

Include them in the YAML front matter:
```yaml
---
title: "Report Title"
metadata-files: ["_fgcz-report.yml"]
---
```

There is **no external CSS file** — all styling lives in `fgcz_header_quarto.html` as inline `<style>` blocks, plus `theme: cosmo` from `_fgcz-report.yml`.

---

## Overall Layout

The entire report body is wrapped in a single top-level `::: {.panel-tabset}`.
Each `# H1` heading becomes a tab. Close with `:::` at the very end.

```markdown
::: {.panel-tabset}

# Settings
...

# Results
...

# Session Info
...

:::
```

---

## Nested Tabsets

Within a top-level tab, use `::: {.panel-tabset}` + `## H2` headings for sub-tabs.
Can nest up to 3 levels deep (H1 > H2 > H3).

```markdown
# Quality Control

::: {.panel-tabset}

## Missing Values
...

## PCA
...

:::
```

---

## Collapsible Callouts

Use fenced div callouts for supplementary information. Always set `collapse="true"`.

| Callout type | Use for |
|---|---|
| `.callout-note` | Analysis parameters, method descriptions, thresholds |
| `.callout-tip` | Links (B-Fabric, downloads), overview summaries |
| `.callout-warning` | Caveats, interpretation warnings, technical bias notes |

```markdown
::: {.callout-note collapse="true"}
## Analysis Parameters
Content here (tables, text, etc.)
:::
```

---

## Tables

**All tables use `DT::datatable()`** — never `kable()` or `kableExtra`.

Standard wrapper pattern:
```r
DT::datatable(
  df, caption = "...", filter = "none",
  rownames = FALSE, class = "compact stripe",
  options = list(pageLength = 10, scrollX = TRUE, dom = "tip")
)
```

---

## Figures

- Use Quarto chunk options (`#| fig-width:`, `#| fig-height:`) not knitr `fig.width`
- `_fgcz-report.yml` sets `lightbox: true` globally — all figures get click-to-zoom
- Save plots to `plots/` as PNG + PDF for download:
  ```r
  ggplot2::ggsave(file.path("plots", "name.png"), plot, width = 7, height = 5, dpi = 300)
  ggplot2::ggsave(file.path("plots", "name.pdf"), plot, width = 7, height = 5)
  ```

---

## Programmatic Tab Generation

When the number of tabs is data-driven (e.g. per-contrast, per-category), use
`#| results: asis` + `cat()` to emit tabset markup:

```r
#| echo: false
#| results: asis
cat("::: {.panel-tabset}\n\n")
for (name in names(results)) {
  cat("##", name, "\n\n")
  print(plot_for(name))
  cat("\n\n")
}
cat(":::\n\n")
```

---

## R Code Conventions

- `#| include: false` for setup/data-loading chunks (no output)
- `#| echo: false` for output-only chunks (plots, tables)
- Use `library()` calls in setup (not explicit namespacing) — this is the report convention, not package code
- Conditional evaluation via `#| eval: !expr "condition"`

---

## Defaults from `_fgcz-report.yml`

These are set globally and should not be overridden per-chunk unless necessary:

| Setting | Value |
|---|---|
| `embed-resources` | `true` (self-contained HTML) |
| `theme` | `cosmo` |
| `code-fold` | `true` |
| `lightbox` | `true` |
| `fig-dpi` | 300 |
| `execute.echo` | `false` |
| `execute.warning` | `false` |
| `knitr.opts_chunk.out.width` | `"30%"` |
| `knitr.opts_chunk.fig.retina` | 2 |
| Grid: body-width | 1800px |

---

## Report Skeleton

Minimal template for a new FGCZ report:

````markdown
---
title: "My Report"
metadata-files: ["_fgcz-report.yml"]
---

```{r setup}
#| include: false
library(DT)
library(ggplot2)
# ... load data ...
```

::: {.panel-tabset}

# Settings

::: {.callout-note collapse="true"}
## Parameters
```{r}
DT::datatable(params_df)
```
:::

# Results

::: {.panel-tabset}

## Overview
...

## Details
...

:::

# Session Info

```{r}
sessionInfo()
```

:::
````