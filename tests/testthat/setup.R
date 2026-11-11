# Run-wide test setup (testthat 3e runs this once before the test files).
#
# Plots that query or open a device implicitly (par("ps"), cowplot legends,
# recordPlot) would otherwise open the default file device and leave
# Rplots.pdf in tests/testthat/ after every run. A null PDF device takes those
# draws instead and is closed when the run ends.
grDevices::pdf(NULL)
withr::defer(grDevices::dev.off(), teardown_env())
