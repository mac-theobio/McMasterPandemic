.onLoad <- function(lib, pkg) {
    flex_version <- readLines(
        system.file("tmb/recommended_spec_version",
            package = "McMasterPandemic"
        )
    )[1]
    options(

        # -- spec version settings ---------------------------------------------

        # default spec version
        MP_flex_spec_version = flex_version,

        # url containing spec version descriptions
        MP_flex_spec_doc_site = "https://canmod.net/misc/flex_specs",

        # default object file for c++ implementation of the spec
        #   - the default means that src/McMasterPandemic.src will be used
        #   - a common alternative is "macpan", which will correspond to
        #     different spec versions in inst/tmb
        MP_flex_spec_dll = "McMasterPandemic",

        # https://stackoverflow.com/questions/8396577/check-if-character-value-is-a-valid-r-object-name
        MP_name_search_regex = "((([[:alpha:]]|[.][._[:alpha:]])[._[:alnum:]]*)|[.])",

        # -- default behaviour of classic macpan engine ------------------------

        # ??
        MP_badsum_action = "warning",

        # ??
        MP_badsum_tol = 1e-12,

        # default number of steps to take in the iterative eigenvector solver
        MP_rexp_steps_default = 100,

        # -- warnings associated with flexmodel structure ----------------------

        # warn if there are multiple rates for the same state transition
        MP_warn_repeated_rates = FALSE,

        # warn if there is no outflow in the model (rarely what users want)
        # -- setting to FALSE because of MP_auto_outflow
        MP_warn_no_outflow = FALSE,

        # warn if there are time-variation breakpoints outside of the
        # simulation bounds
        MP_warn_bad_breaks = TRUE,

        # -- flexmodel defaults ------------------------------------------------
        # -- see ?flexmodel for description of these arguments -----------------
        MP_default_do_hazard = FALSE,
        MP_default_do_hazard_lin = FALSE,
        MP_default_do_approx_hazard = FALSE,
        MP_default_do_approx_hazard_lin = FALSE,
        MP_default_do_make_state = FALSE,
        MP_default_do_sim_constraint = FALSE,
        MP_default_sim_lower_bound = 1e-12,

        # -- tmb_fun behaviour -------------------------------------------------

        # if the user does not specify any outflows, automatically add
        # an outflow for every inflow
        MP_auto_outflow = TRUE,

        # automatically update the tmb indices when (re)generating the
        # tmb ad fun
        MP_auto_tmb_index_update = TRUE,

        # do state condensation with c++ code as opposed to r
        MP_condense_cpp = TRUE,

        # -- control how comparable r and tmb engines are ----------------------
        # -- see r_tmb_comparable ----------------------------------------------
        MP_use_state_rounding = TRUE,         # FALSE ~ comparable
        MP_vax_make_state_with_hazard = TRUE, # FALSE ~ comparable
        MP_tmb_models_match_r = FALSE,        # TRUE  ~ comparable
        MP_force_dgTMatrix = FALSE            # TRUE  ~ comparable
    )
}

.onUnload <- function(libpath) {
    library.dynam.unload("McMasterPandemic", libpath)
}
