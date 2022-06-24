remotes::install_github("canmod/macpan2")
library(McMasterPandemic)
library(macpan2)
parser = make_expr_parser('parser', finalizer_index)
model = (make_hello_world_model()
  %>% update_simulation_bounds('2000-01-01', '2000-01-02')
)
valid_funcs = McMasterPandemic:::nlist(`+`, `-`, `*`, `/`, `^`, `(`)
valid_vars = as.list(c(model$state, model$params))

test_expr = ~ (100 * beta * I / N) + S/100

expr_parsed = parser(test_expr)
expr_parsed$parse_table
names(expr_parsed$valid_funcs)
names(expr_parsed$valid_vars)
unlist(expr_parsed$valid_literals)

eval(test_expr[[2]], expr_parsed$valid_vars)
with(simulation_history(model)[1, ], 100 * S_to_I + S/100)

args = tmb_fun_args(model)
args$data$xx = expr_parsed$parse_table$x
args$data$nn = expr_parsed$parse_table$n
args$data$ii = expr_parsed$parse_table$i
args$data$vv = unname(unlist(expr_parsed$valid_vars))
args$data$ll = unname(unlist(expr_parsed$valid_literals))
fn = do.call(MakeADFun, args)
fn$fn() # this should equal the following two expressions

eval(test_expr[[2]], expr_parsed$valid_vars)
with(simulation_history(model)[1, ], 100 * S_to_I + S/100)
