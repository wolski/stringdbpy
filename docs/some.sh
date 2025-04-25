# convert you ipynb to qmd
quarto convert docs/python_notebooks/VisualizeNetworks.ipynb  -o docs/VisualizeNetworks.qmd

# in docs directory
quarto render . 
quarto preview .
