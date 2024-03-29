#' @importFrom dplyr mutate select left_join
#' @importFrom brms brm mvbf bf prior
#' @importFrom ergm ergm control.ergm
#' @importFrom purrr map_dfr
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stats as.formula coef
#' @importFrom utils head tail
#' @importFrom rlang enquo quo_name .data
#' @importFrom magrittr %>%
#' @importFrom network %v%