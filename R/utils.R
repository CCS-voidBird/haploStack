# =============================================================================
# utils.R — small helpers with no owner module
# =============================================================================

# Null-coalescing helper (not exported)
`%||%` <- function(a, b) if (is.null(a)) b else a


# Reorder a data.table in place so its identification columns lead. `priority`
# lists id columns most-significant first; those present are pulled to the front
# in that order and the remaining columns follow in their current relative order.
# Returns the table invisibly (modified in place via setcolorder).
# @keywords internal
.lead_identifier_cols <- function(dt, priority = .HS_ID_COLS) {
  nm    <- names(dt)
  lead  <- priority[priority %in% nm]
  rest  <- setdiff(nm, lead)
  data.table::setcolorder(dt, c(lead, rest))
  invisible(dt)
}

