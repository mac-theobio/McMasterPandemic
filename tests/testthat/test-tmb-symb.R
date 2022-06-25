remotes::install_github("canmod/macpan2")
library(McMasterPandemic)
library(macpan2)
set_spec_version('0.3.0', 'inst/tmb')
parser = make_expr_parser('parser', finalizer_index)
model = (make_hello_world_model()
  %>% update_simulation_bounds('2000-01-01', '2000-01-02')
  %>% update_observed(data.frame(date = "2000-01-02", var = "S", value = 20000))
  %>% add_opt_params(beta ~ flat(0.15))
)
valid_funcs = McMasterPandemic:::nlist(`+`, `-`, `*`, `/`, `^`, `(`)
valid_vars = list(
  beta = model$params[["beta"]],
  N = model$params[["N"]],
  I = model$state[["I"]],
  S = model$state[["S"]]
)

test_expr = ~ (100 * beta * I / N) + S/100

expr_parsed = parser(test_expr)
expr_parsed$parse_table
names(expr_parsed$valid_funcs)
names(expr_parsed$valid_vars)
unlist(expr_parsed$valid_literals)

args = tmb_fun_args(model)
args$data$parse_table_x = as.integer(expr_parsed$parse_table$x)
args$data$parse_table_n = as.integer(expr_parsed$parse_table$n)
args$data$parse_table_i = as.integer(expr_parsed$parse_table$i)
args$parameters$valid_vars = unname(unlist(expr_parsed$valid_vars))
args$data$valid_literals = unname(unlist(expr_parsed$valid_literals))
params_map = factor(
  rep(NA, length(model$params)),
  levels = c()
)
tv_mult_map = factor(
  rep(NA, nrow(model$timevar$piece_wise$schedule)),
  levels = c()
)
valid_vars_map = factor(
  1:4,
  levels = 1:4
)
args$map = list(
  params = params_map,
  tv_mult = tv_mult_map,
  valid_vars = valid_vars_map
)
fn = do.call(MakeADFun, args)
sims = fn$report()$simulation_history[1, ]

fn$fn() # this should equal the following two expressions
eval(test_expr[[2]], expr_parsed$valid_vars)
sims[1] / 100 + 100 * sims[4]
