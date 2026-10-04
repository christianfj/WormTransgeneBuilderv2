test_that("software_versions reports the app, R and ViennaRNA versions", {
  v <- software_versions()
  expect_named(v, c("app", "R", "ViennaRNA"))
  expect_identical(v[["app"]], as.character(packageVersion("WormTransgeneBuilder")))
  expect_identical(v[["R"]], paste(R.version$major, R.version$minor, sep = "."))
})

test_that("a missing or broken RNAfold is reported as not installed", {
  expect_identical(software_versions(rnafold = "")[["ViennaRNA"]], "not installed")
  expect_identical(software_versions(rnafold = "/no/such/RNAfold")[["ViennaRNA"]], "not installed")
})

test_that("the ViennaRNA version comes from RNAfold itself", {
  skip_if(!nzchar(find_rnafold()), "RNAfold is not installed")
  expect_match(software_versions()[["ViennaRNA"]], "^[0-9]+\\.[0-9]+\\.[0-9]+")
})

test_that("the About tab shows the versions", {
  shiny::testServer(app_server, {
    session$setInputs(checkboxRibo = FALSE)
    expect_match(output$software_versions, "^WormTransgeneBuilder .*, R [0-9.]+, ViennaRNA ")
  })
})
