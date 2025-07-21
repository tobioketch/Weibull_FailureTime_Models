# Directory
THIS_PATH__ <- base::file.path("2--run_codes", "NOT_RUN_RCODES")

# Load packages and functions
source(file = here::here(THIS_PATH__, "MethodsAndPackages.R")) 
# just in case parallel there's a memory leak due to open parallel processing sockets 
future::plan(future::sequential)

# Load sampling models
source(file = here::here(THIS_PATH__, "StanSamplingModels","PriorSamplingSchemes.R") )


