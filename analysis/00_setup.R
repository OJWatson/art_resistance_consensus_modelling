# Setup For correctly using python in this R reproducible compendium package

# Ensure all R packages are installed that are needed
devtools::install_dev_deps()
devtools::load_all(".")

# Installing any python packages not available as default by reticulate
reticulate::py_install("seaborn")
