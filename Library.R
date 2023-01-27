load_packages <- function(packages, quietly=TRUE){
  for(package in packages){
    installed = 
      if(quietly)
        suppressWarnings(suppressMessages(suppressPackageStartupMessages(
          require(package, character.only = TRUE, quietly = TRUE)
        )))
    else
      require(package, character.only = TRUE)
    
    if(!installed){
      install.packages(package,repos = "https://cloud.r-project.org")
      if(quietly)
        suppressWarnings(suppressMessages(suppressPackageStartupMessages(
          library(package, character.only = TRUE, quietly = TRUE)
        )))
      else
        library(package, character.only = TRUE)
    }
  }
  invisible(NULL)     
}

dir_create <- function(paths){
  # only creates the directory if it doesnt exist
  # defaults to recursive
  for(path in paths)
    if(!file.exists(path)) dir.create(path, recursive=TRUE)
  
  invisible(NULL)
}

