# =============================================================================
# test-090-agent-guide.R — the usage guide shipped in inst/skills (0.9.0)
#
#   G1  hs_agent_guide() finds and prints SKILL.md (frontmatter intact)
#   G2  every reference topic prints; an unknown topic is an input error
#   G3  hs_install_skill() copies the skill, stamps the version, refuses to
#       overwrite unless asked
#   G4  the guide stays true: every function it calls exists (a haploStack
#       export, or base / recommended R), every topic it names exists, ASCII
#   G5  the startup line names hs_agent_guide()
# =============================================================================

.g090_pages <- function() {
  d <- system.file("skills", "haplostack", package = "haploStack")
  c(file.path(d, "SKILL.md"),
    list.files(file.path(d, "references"), pattern = "[.]md$", full.names = TRUE))
}

test_that("G1: hs_agent_guide() prints SKILL.md and returns its path", {
  out <- capture.output(f <- hs_agent_guide())
  expect_true(file.exists(f))
  expect_identical(basename(f), "SKILL.md")
  expect_identical(out[1], "---")                       # frontmatter first
  expect_true(any(out == "name: haplostack"))
  expect_true(any(grepl("^description: ", out)))
  expect_identical(capture.output(g <- hs_agent_guide(print = FALSE)), character(0))
  expect_identical(g, f)
})

test_that("G2: every reference topic prints; unknown topics are refused", {
  topics <- haploStack:::.hs_skill_topics()
  expect_true(all(c("blocks", "genometa", "simulation", "troubleshooting") %in% topics))
  for (t in topics) {
    out <- capture.output(f <- hs_agent_guide(t))
    expect_identical(basename(f), paste0(t, ".md"))
    expect_gt(length(out), 5L)
  }
  e <- tryCatch(hs_agent_guide("nope"), error = identity)
  expect_s3_class(e, "hs_error_input")
  expect_match(conditionMessage(e), "\"blocks\"", fixed = TRUE)
  expect_s3_class(tryCatch(hs_agent_guide(c("blocks", "genometa")), error = identity),
                  "hs_error_input")
})

test_that("G3: hs_install_skill() copies, stamps and guards the target", {
  dest <- tempfile("skills_")
  on.exit(unlink(dest, recursive = TRUE), add = TRUE)
  expect_message(d <- hs_install_skill(dest), "installed in", class = "hs_message")
  expect_identical(normalizePath(d), normalizePath(file.path(dest, "haplostack")))
  src <- system.file("skills", "haplostack", package = "haploStack")
  rel <- function(x) sort(list.files(x, recursive = TRUE))
  expect_identical(rel(d), rel(src))
  txt <- readLines(file.path(d, "SKILL.md"))
  stamp <- paste0("installed from haploStack ", utils::packageVersion("haploStack"))
  expect_true(any(grepl(stamp, txt, fixed = TRUE)))
  expect_false(any(grepl("<!-- haplostack-version -->", txt, fixed = TRUE)))
  expect_identical(txt[1], "---")                       # frontmatter still first
  # the shipped copy keeps the placeholder
  expect_true(any(grepl("<!-- haplostack-version -->",
                        readLines(file.path(src, "SKILL.md")), fixed = TRUE)))
  # a second install is refused, overwrite = TRUE replaces
  writeLines("stale", file.path(d, "stale.txt"))
  expect_s3_class(tryCatch(hs_install_skill(dest), error = identity), "hs_error_input")
  expect_true(file.exists(file.path(d, "stale.txt")))
  hs_quiet(hs_install_skill(dest, overwrite = TRUE))
  expect_false(file.exists(file.path(d, "stale.txt")))
  expect_identical(rel(d), rel(src))
  # 'path' is required and must be one name
  expect_s3_class(tryCatch(hs_install_skill(), error = identity), "hs_error_input")
  expect_s3_class(tryCatch(hs_install_skill(c("a", "b")), error = identity),
                  "hs_error_input")
})

test_that("G4: every function the guide calls exists; every topic it names exists", {
  pages <- .g090_pages()
  expect_gte(length(pages), 5L)
  txt <- unlist(lapply(pages, readLines, warn = FALSE))
  # ASCII only: the guide is printed to consoles of any locale
  expect_false(any(grepl("[^\\x01-\\x7F]", txt, perl = TRUE)),
               label = paste("non-ASCII lines:",
                             paste(head(grep("[^\\x01-\\x7F]", txt, perl = TRUE,
                                             value = TRUE), 3), collapse = " | ")))
  # bare calls name(  — not pkg::name(, not $name(, not .name(
  m <- regmatches(txt, gregexpr("(?<![[:alnum:]_.:$])[A-Za-z][A-Za-z0-9._]*(?=\\()",
                                txt, perl = TRUE))
  calls <- sort(unique(unlist(m)))
  exports <- getNamespaceExports("haploStack")
  other <- c("base", "stats", "utils", "graphics", "grDevices", "methods",
             "data.table")
  known <- function(f) f %in% exports ||
    any(vapply(other, function(p) exists(f, envir = asNamespace(p), inherits = FALSE),
               logical(1)))
  # names of other packages' functions the guide mentions without a prefix
  foreign <- c("st.def.hblocks")
  unknown <- setdiff(calls[!vapply(calls, known, logical(1))], foreign)
  expect_identical(unknown, character(0),
                   label = paste("functions the guide names that do not exist:",
                                 paste(unknown, collapse = ", ")))
  # topics named as hs_agent_guide("x") exist
  named <- unique(unlist(regmatches(txt, gregexpr("(?<=hs_agent_guide\\(\")[a-z]+",
                                                  txt, perl = TRUE))))
  expect_true(length(named) >= 4L)
  expect_true(all(named %in% haploStack:::.hs_skill_topics()),
              label = paste("missing topics:",
                            paste(setdiff(named, haploStack:::.hs_skill_topics()),
                                  collapse = ", ")))
})

test_that("G5: the startup line names hs_agent_guide()", {
  expect_message(haploStack:::.onAttach("", "haploStack"), "hs_agent_guide()",
                 fixed = TRUE)
})
