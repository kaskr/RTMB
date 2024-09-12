if (requireNamespace("tinytest", quietly=TRUE)) {
    tinytest::test_package("RTMB", ncpu=getOption("Ncpus", 1))
}
