# Starts the app inside the container and first prints the versions it runs with.
v <- WormTransgeneBuilder:::software_versions()
cat(sprintf("WormTransgeneBuilder %s | R %s | ViennaRNA %s\n", v[["app"]], v[["R"]], v[["ViennaRNA"]]))
options(shiny.port = 3838, shiny.host = "0.0.0.0")
WormTransgeneBuilder::run_app(options = list(launch.browser = FALSE))
