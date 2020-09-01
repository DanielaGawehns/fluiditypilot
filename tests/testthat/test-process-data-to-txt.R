context("test-process-data")

data <- utils::read.csv(system.file("testdata",
                                    "PilotBeagle.csv",
                                    package = "fluiditypilot"))



test_that("multiplication works", { #test desc
  expect_equal(2 * 2, 4)            #test expectation
})
