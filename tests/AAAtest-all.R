if (require("testthat") & require("here")) {
    print("AAA working directory: ")
    print(getwd())
    print("AAA here: ")
    print(here())
    pkg <- "McMasterPandemic"
    require(pkg, character.only = TRUE)
    print(sessionInfo())
    test_check(pkg, reporter = "summary")
    print(warnings()) ## TODO? catch most of these by expect_warning(..)
} else {
    warnings(
        "Package 'testthat' not available, cannot run unit tests for package",
        sQuote(pkg)
    )
}
