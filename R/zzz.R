.onLoad <- function(lib, pkg) {
    flex_version <- readLines(
        system.file("tmb/recommended_spec_version",
            package = "McMasterPandemic"
        )
    )[1]
    options(

        # spec version settings
        MP_flex_spec_version = flex_version,
        MP_flex_spec_doc_site = "https://canmod.net/misc/flex_specs",
        MP_flex_spec_dll = "McMasterPandemic",

        # https://stackoverflow.com/questions/8396577/check-if-character-value-is-a-valid-r-object-name
        MP_name_search_regex = "((([[:alpha:]]|[.][._[:alpha:]])[._[:alnum:]]*)|[.])",

        # default behaviour of classic macpan engine
        MP_badsum_action = "warning",
        MP_badsum_tol = 1e-12,
        MP_rexp_steps_default = 100,

        # warnings associated with flex model structure
        MP_warn_repeated_rates = FALSE,
        MP_warn_no_outflow = TRUE,
        MP_warn_bad_breaks = TRUE,

        # init_model defaults
        MP_default_do_hazard = TRUE,
        MP_default_do_make_state = TRUE,

        # tmb_fun behaviour
        MP_force_full_outflow = FALSE,
        MP_auto_tmb_index_update = TRUE,

        # control how comparable r and tmb engines are
        # - see r_tmb_comparable
        MP_use_state_rounding = TRUE,         # FALSE ~ comparable
        MP_vax_make_state_with_hazard = TRUE, # FALSE ~ comparable
        MP_tmb_models_match_r = FALSE,         # TRUE ~ comparable
        MP_force_dgTMatrix = FALSE            # TRUE ~ comparable
    )
}

.onUnload <- function(libpath) {
    library.dynam.unload("McMasterPandemic", libpath)
}
