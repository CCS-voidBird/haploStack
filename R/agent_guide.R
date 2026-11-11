# =============================================================================
# agent_guide.R — the usage guide shipped inside the package (0.9.0)
#
# inst/skills/haplostack/ is an agent skill (SKILL.md + references/*.md) that
# installs with the package, so the guide always matches the installed build.
# Coding agents do not search the R library for skills; they find this one
# through the startup line (.onAttach, R/zzz.R), the package help page and the
# README, all of which name hs_agent_guide(). hs_install_skill() copies the
# skill into an agent's skills folder on request — never on install or load.
# test-090-agent-guide.R keeps every function the guide names an export.
# =============================================================================

.hs_skill_dir <- function() {
  d <- system.file("skills", "haplostack", package = "haploStack")
  if (!nzchar(d) || !file.exists(file.path(d, "SKILL.md")))
    .hs_abort("haploStack's usage guide (inst/skills/haplostack) is missing ",
              "from this installation; reinstall the package.",
              class = "hs_error_internal")
  d
}

# Topic names = the reference file names without ".md".
.hs_skill_topics <- function(dir = .hs_skill_dir()) {
  sub("[.]md$", "", .hs_sort(list.files(file.path(dir, "references"), pattern = "[.]md$")))
}

#' Print haploStack's usage guide (for people and coding agents)
#'
#' haploStack ships a usage guide inside the package: an agent skill
#' (\code{SKILL.md}) with reference pages on the block map, the
#' \code{GenoMeta} container, the stacking simulation and troubleshooting. It
#' is installed with the package, so it describes exactly the build in use.
#' \code{hs_agent_guide()} prints it; an AI coding agent that meets haploStack
#' for the first time should run it before anything else.
#'
#' @param topic \code{NULL} (default) for the main guide, or the name of a
#'   reference page: \code{"blocks"}, \code{"genometa"}, \code{"simulation"}
#'   or \code{"troubleshooting"}.
#' @param print Logical; write the page to the console (default \code{TRUE}).
#'   \code{FALSE} only returns the path.
#' @return The path of the page, invisibly.
#' @seealso \code{\link{hs_install_skill}} to copy the guide into an agent's
#'   skills folder.
#' @examples
#' hs_agent_guide(print = FALSE)            # where the guide lives
#' \dontrun{
#' hs_agent_guide()                         # the main guide
#' hs_agent_guide("blocks")                 # one reference page
#' }
#' @export
hs_agent_guide <- function(topic = NULL, print = TRUE) {
  dir <- .hs_skill_dir()
  if (is.null(topic)) {
    f <- file.path(dir, "SKILL.md")
  } else {
    topics <- .hs_skill_topics(dir)
    if (!is.character(topic) || length(topic) != 1L || !topic %in% topics)
      .hs_abort("hs_agent_guide: 'topic' must be NULL or one of ",
                paste0("\"", topics, "\"", collapse = ", "), ".")
    f <- file.path(dir, "references", paste0(topic, ".md"))
  }
  if (isTRUE(print)) writeLines(readLines(f, warn = FALSE, encoding = "UTF-8"))
  invisible(f)
}

#' Copy haploStack's usage guide into an agent's skills folder
#'
#' Copies the agent skill that ships with the package
#' (\code{inst/skills/haplostack}, see \code{\link{hs_agent_guide}}) to
#' \code{file.path(path, "haplostack")}, so that an agent which loads skills
#' from that folder finds it in later sessions without being told. Nothing is
#' written unless this function is called. Typical folders: a project's
#' \code{.claude/skills} (Claude Code) or \code{.agents/skills}, or a user-level
#' folder such as \code{~/.claude/skills}.
#'
#' The copy records the haploStack version it came from in a comment line of
#' its \code{SKILL.md}. After upgrading haploStack, refresh it with
#' \code{overwrite = TRUE}.
#'
#' @param path Folder that holds the agent's skills. Created if absent.
#' @param overwrite Logical; replace an existing \code{haplostack} skill in
#'   \code{path} (default \code{FALSE}, which refuses).
#' @return The installed skill folder, invisibly.
#' @examples
#' dest <- file.path(tempdir(), "skills")
#' hs_install_skill(dest)
#' list.files(file.path(dest, "haplostack"), recursive = TRUE)
#' @export
hs_install_skill <- function(path, overwrite = FALSE) {
  if (missing(path) || !is.character(path) || length(path) != 1L ||
      is.na(path) || !nzchar(path))
    .hs_abort("hs_install_skill: 'path' must be one folder name, e.g. ",
              "\".claude/skills\".")
  src  <- .hs_skill_dir()
  dest <- file.path(path, "haplostack")
  if (dir.exists(dest)) {
    if (!isTRUE(overwrite))
      .hs_abort("hs_install_skill: '", dest, "' already exists. Pass ",
                "overwrite = TRUE to replace it.")
    unlink(dest, recursive = TRUE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  # copy the folder's contents, not the folder, so the target name is fixed
  dir.create(dest)
  ok <- file.copy(list.files(src, full.names = TRUE), dest, recursive = TRUE)
  if (!all(ok))
    .hs_abort("hs_install_skill: could not copy the guide into '", dest, "'.")
  skill <- file.path(dest, "SKILL.md")
  txt <- readLines(skill, warn = FALSE, encoding = "UTF-8")
  version <- as.character(utils::packageVersion("haploStack"))
  txt <- sub("<!-- haplostack-version -->",
             paste0("<!-- installed from haploStack ", version, " on ",
                    format(Sys.Date()), "; refresh with ",
                    "haploStack::hs_install_skill(path, overwrite = TRUE) ",
                    "after an upgrade -->"),
             txt, fixed = TRUE)
  con <- file(skill, open = "w", encoding = "UTF-8")
  writeLines(txt, con)
  close(con)
  .hs_inform("haplostack skill (haploStack ", version, ") installed in ", dest)
  invisible(dest)
}
