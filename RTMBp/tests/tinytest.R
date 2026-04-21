if (requireNamespace("tinytest", quietly=TRUE)) {
    tinytest::test_package("RTMBp", ncpu=getOption("Ncpus", 1))
}
