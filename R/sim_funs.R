##' Compute elementwise multiplication of a vector by each column in a matrix
##' This is implemented for performance reasons, as it is much faster than sweep.
##' @param M matrix
##' @param v vector
##' @export
col_multiply <- function(M, v) {
    new <- M * v ## t(t(M) * rep(v, rep.int(nrow(M), length(v))))
    ##    old <- sweep(M, v, MARGIN=1, FUN="*")
    ##    print(all(new==old))
    return(new)
}

##' construct Jacobian matrix for ICU model
##' (not quite complete: doesn't include flows to R)
## FIXME: derive from make_ratemat
##' @param state state vector (named)
##' @param params parameter vector
##' @export
##' @examples
##' params <- read_params("ICU1.csv")
##' state <- make_state(params[["N"]],E0=params[["E0"]], use_eigvec=FALSE)
##' ## state[c("E","Ia","Ip","Im","Is")] <- 1
##' state[["E"]] <- 1
##' J <- make_jac(params,state)
##' J["S","S"]
##' Jr <- J[1:6,1:6]
##' round(Jr,3)
##' eigen(Jr)$values
##' make_jac(params)
make_jac <- function(params, state = NULL) {
    ## circumvent test code analyzers ... problematic ...
    S <- E <- Ia <- Ip <- Im <- Is <- H <- NULL
    H2 <- ICUs <- ICUd <- D <- R <- beta0 <- Ca <- Cp <- NULL
    Cm <- Cs <- alpha <- sigma <- gamma_a <- gamma_m <- gamma_s <- gamma_p <- NULL
    rho <- delta <- mu <- N <- E0 <- iso_m <- iso_s <- phi1 <- NULL
    phi2 <- psi1 <- psi2 <- psi3 <- c_prop <- c_delaymean <- c_delayCV <- NULL
    ##
    if (is.null(state)) {
        state <- make_state(
            N = params[["N"]], E0 = 1e-3,
            use_eigvec = FALSE
        )
    }
    np <- length(params)
    ns <- length(state)
    ## make state and param names locally available (similar to with())
    P <- c(as.list(state), as.list(params))
    unpack(P) ## extract variables
    ## blank matrix
    M <- matrix(0,
        nrow = ns, ncol = ns,
        dimnames = list(from = names(state), to = names(state))
    )
    Ivec <- c(Ia, Ip, Im, Is)
    Iwt <- beta0 / N * c(Ia = Ca, Ip = Cp, Im = (1 - iso_m) * Cm, Is = (1 - iso_s) * Cs)
    Ivars <- c("Ia", "Ip", "Im", "Is")
    M["S", "S"] <- -sum(Ivec * Iwt)
    M["S", Ivars] <- -S * Iwt[Ivars]
    M["E", c("S", Ivars)] <- -M["S", c("S", Ivars)]
    M["E", "E"] <- -sigma
    M["Ia", "E"] <- alpha * sigma
    M["Ia", "Ia"] <- -gamma_a
    M["Ip", "E"] <- (1 - alpha) * sigma
    M["Ip", "Ip"] <- -gamma_p
    M["Im", "Ip"] <- mu * gamma_p
    M["Im", "Im"] <- -gamma_m
    M["Is", "Ip"] <- (1 - mu) * gamma_p
    M["Is", "Is"] <- -gamma_s
    M["H", "Is"] <- phi1 * gamma_s
    M["H", "H"] <- -rho
    M["ICUs", "Is"] <- (1 - phi1) * (1 - phi2) * gamma_s
    M["ICUs", "ICUs"] <- -psi1
    M["H2", "ICUs"] <- psi1
    M["H2", "H2"] <- -psi3
    M["ICUd", "Is"] <- (1 - phi1) * phi2 * gamma_s
    M["ICUd", "ICUd"] <- -psi2
    M["D", "ICUd"] <- psi2
    M["R", "Ia"] <- gamma_a
    M["R", "Im"] <- gamma_m
    M["R", "H"] <- rho
    M["R", "H2"] <- psi3
    return(M)
}

##' construct vector of transmission multipliers
##' @param state state vector
##' @param params parameter vector
##' @param full include non-infectious compartments (with transmission of 0) as well as infectious compartments?
##' @export
## QUESTION: is the main testify argument to this function used?
make_betavec <- function(state, params, full = TRUE) {
  Icats <- c("Ia", "Ip", "Im", "Is")
  testcats <- c("_u", "_p", "_n", "_t")
  ## NB meaning of iso_* has switched from Stanford model
  ## beta_vec0 is the vector of transmission parameters that apply to infectious categories only
  beta_vec0 <- with(
    as.list(params),
    beta0 / N * c(Ca, Cp, (1 - iso_m) * Cm, (1 - iso_s) * Cs)
  )
  names(beta_vec0) <- Icats
  if (has_age(params)) {
    Cmat <- params$Cmat
    a_names <- rownames(Cmat)
    new_names <- expand_names(Icats, a_names)
    beta_vec0 <- t(kronecker(Cmat, matrix(beta_vec0)))
    dimnames(beta_vec0) <- list(a_names, new_names)
    beta_vec0 <- Matrix(beta_vec0)
  }
  ## assume that any matching values will be of the form "^%s_" where %s is something in Icats
  ## lapply(Icats, function(x) grep(sprintf("^%s_"), names(state))
  ## FIXME: we should be doing this by name, not assuming that all infectious compartments are expanded
  ##  into exactly 4 subcompartments, in order (but this should work for now??)
  if (has_testing(state = state)) { ## testified!
    if (has_age(params)) stop("can't combine age and testing yet")
    beta_vec0 <- rep(beta_vec0, each = length(testcats))
    names(beta_vec0) <- unlist(lapply(Icats, function(x) paste0(x, testcats)))
    ## FIXME: also adjust _n, _p components?
    pos_vals <- grep("_t$", names(beta_vec0))
    beta_vec0[pos_vals] <- beta_vec0[pos_vals] * (1 - params[["iso_t"]])
  }
  if (!full) {
    return(beta_vec0)
  }
  ## By default, make a vector of zeroes for all the states,
  ## then fill in infectious ones
  if (!has_age(params)) {
    beta_vec <- setNames(numeric(length(state)), names(state))
    beta_vec[names(beta_vec0)] <- beta_vec0
  } else {
    beta_vec <- matrix(0,
                       nrow = nrow(beta_vec0), ncol = length(state),
                       dimnames = list(rownames(beta_vec0), names(state))
    )
    beta_vec[rownames(beta_vec0), colnames(beta_vec0)] <- matrix(beta_vec0)
  }
  return(beta_vec)
}

calc_variant_adjustment <- function(params,
                                    stratum = c("unvax", "dose1", "dose2")) {
    ## used for the vaxified model
    stratum <- match.arg(stratum)

    ## calculate adjustments to beta0 based on increased transmissibility and/or vaccine efficacy
    if (stratum == "unvax") {
        dominant_efficacy <- variant_efficacy <- 0
    } else {
        dominant_efficacy <- params[[paste0("vax_efficacy_", stratum)]]
        variant_efficacy <- params[[paste0("variant_vax_efficacy_", stratum)]]
    }

    adjustment <- (1 - dominant_efficacy) * (1 - params[["variant_prop"]]) + (1 - variant_efficacy) * params[["variant_advantage"]] * params[["variant_prop"]]

    return(adjustment)
}

##' construct vector of transmission multipliers
##' @param state state vector
##' @param params parameter vector
##' @param full include non-infectious compartments (with transmission of 0) as well as infectious compartments?
##' @importFrom fastmatrix kronecker.prod
##' @export
## QUESTION: is the main testify argument to this function used?
make_beta <- function(state, params, full = TRUE) {
    Icats <- c("Ia", "Ip", "Im", "Is")
    testcats <- c("_u", "_p", "_n", "_t")
    ## NB meaning of iso_* has switched from Stanford model
    ## beta_vec0 is the vector of transmission parameters
    ## that apply to infectious categories only
    ##
    ## (delaying normalization by population size here in
    ## case the age-structured model is being used, where we
    ## normalize by the size of the age-specific susceptible
    ## population involved instead)
    Icat_prop_vec <- with(
        as.list(params),
        c(Ca, Cp, (1 - iso_m) * Cm, (1 - iso_s) * Cs)
    )
    names(Icat_prop_vec) <- Icats

    ## initialize beta_0, which includes combines all parameters for FOI term
    ## (except for state values(
    ## without age, this is beta0*Icat_prop_vec/N
    ## with age, we additionally have contact probabilities stored in pmat
    ## (and the fact that N and maybe beta0 are age-specific)
    if (has_age(params)) {

        ## check that all components for ageified transmission multipliers are
        ## present in params and in the right format
        if (!is.list(params)) stop("must expand parameters to include age componets first (use expand_params_age())")

        ## pmat checks
        if (is.null(params$pmat)) stop("must specify params$pmat component")
        ## check that pmat rows sum to 1
        if (!isTRUE(all.equal(unname(rowSums(params$pmat)), rep(1, nrow(params$pmat))))) stop("each pmat row must sum to 1 (it should be a probability distribution)")

        ## Nvec check
        if (length(params$N) != nrow(params$pmat)) stop("N must be a vector of the same length as the number of age groups specified via pmat.")

        ## beta0 check
        if (!(length(params$beta0 %in% c(1, nrow(params$pmat))))) stop("beta0 must either be a scalar (same beta0 for all ages) or a vector of the same length as the number of age groups specified via pmat.")

        ## incorporate contact matrix and /N_j in beta term, and attach age cats

        ## grab contact matrix (with susceptibles as rows and infectives as
        ## columns) and scale each row by beta0 corresponding to that
        ## susceptible age group (if beta0 is a scalar, scale the whole matrix
        ## with it)
        pmat <- params$beta0 * params$pmat
        ## transpose newly-scaled pmat to enable calculations below
        ##
        ## calculate c_{ij}/N_j to incorporate 1/N_j from
        ## I_j/N_j in force of infection (for mat/vec, R will
        ## divide the entire first row of mat by the first
        ## element of vec, etc.)
        pmat <- t(pmat) / params$N
        a_names <- rownames(pmat)
        new_names <- expand_names(Icats, a_names)
        ## transpose back so rows represent susceptibles and
        ## columns represent infectives again
        beta_0 <- Matrix::t(kronecker(pmat, Matrix::Matrix(Icat_prop_vec)))
        dimnames(beta_0) <- list(a_names, new_names)
    } else {
        ## without age structure, multiply by single beta0 and
        ## normalize by total population N for I/N term in force of infection
        beta_0 <- with(as.list(params), beta0 * Icat_prop_vec / N)
    }

    ## handle vaccination, if present
    if (has_vax(params)) {
        if (!has_vax(state)) stop("if params are vaxified, state also needs to be vaxified")
        ## get vax categories
        vax_cat <- get_vax(params)
        model_type <- get_vax_model_type(vax_cat)

        ## save original beta_0 names
        ## of beta_0 is just a vector (base case), just get colnames
        if (is.null(dimnames(beta_0))) {
            original_row_names <- NULL
            original_col_names <- names(beta_0)
        } else {
            ## otherwise, get both row and colnames
            original_row_names <- rownames(beta_0)
            original_col_names <- colnames(beta_0)
        }

        ## initialize vaccine transmission reduction matrix
        ## for unvax and vaxdose categories, assume no changes to transmission
        ## for vaxprotect categories, assume reduction to transmission equivalent
        ## to vaccine efficacy (for the specified dose)
        vax_trans_red <- matrix(1,
            nrow = length(vax_cat),
            ncol = length(vax_cat)
        )
        rownames(vax_trans_red) <- vax_cat

        ## incorporate vaccine efficacy in different ways, depending on whether or not there is a variant
        if (!do_variant(params)) {
            vax_trans_red[grepl(vax_cat[3], rownames(vax_trans_red)), ] <- rep(1 - params[["vax_efficacy_dose1"]], length(vax_cat))

            if (model_type == "twodose") {
                ## carry over efficacy from one dose to vaxdose2 (received second dose, but not yet protected)
                vax_trans_red[grepl(vax_cat[4], rownames(vax_trans_red)), ] <- rep(1 - params[["vax_efficacy_dose1"]], length(vax_cat))

                ## add in efficacy for second dose
                vax_trans_red[grepl(vax_cat[5], rownames(vax_trans_red)), ] <- rep(1 - params[["vax_efficacy_dose2"]], length(vax_cat))
            }
        } else {
            ## if we're including a variant, we also need to take care to do the increased transmissibility adjustment here (for non-vaxified models, that's taken care of below)
            vax_trans_red[grepl(vax_cat[1], rownames(vax_trans_red)), ] <- rep(calc_variant_adjustment(params, "unvax"))
            vax_trans_red[grepl(vax_cat[2], rownames(vax_trans_red)), ] <- rep(calc_variant_adjustment(params, "unvax"))

            vax_trans_red[grepl(vax_cat[3], rownames(vax_trans_red)), ] <- rep(calc_variant_adjustment(params, "dose1"))

            if (model_type == "twodose") {
                vax_trans_red[grepl(vax_cat[4], rownames(vax_trans_red)), ] <- rep(calc_variant_adjustment(params, "dose1"))

                vax_trans_red[grepl(vax_cat[5], rownames(vax_trans_red)), ] <- rep(calc_variant_adjustment(params, "dose2"))
            }
        }

        ## apply vaccine transmission reduction over beta_0 (as computed above,
        ## either with or without age) using the kronecker product trick
        if (class(beta_0) == "numeric") {
            ## need to take the transpose of the kronecker result since
            ## Matrix::Matrix(beta_0) in the kronecker product takes a row vector
            ## and converts it to a column vector
            ##    beta_0_old <- kronecker(Matrix::Matrix(vax_trans_red),
            ##                        Matrix::t(Matrix::Matrix(beta_0)))
            beta_0_fastmatrix <- Matrix::Matrix(fastmatrix::kronecker.prod(vax_trans_red, t(beta_0)))
            ##    beta_0_base <- Matrix::Matrix(vax_trans_red %x% t(beta_0))
            ##    beta_0_rtensor <- Matrix::Matrix(rTensor::kronecker_list(list(vax_trans_red, t(beta_0))))
        } else {
            ## otherwise, beta_0 will already be a Matrix::Matrix (class dgeMatrix) and we don't need any transposes
            ##      beta_0_old <- kronecker(Matrix::Matrix(vax_trans_red), Matrix::Matrix(beta_0))

            beta_0_fastmatrix <- Matrix::Matrix(fastmatrix::kronecker.prod(vax_trans_red, as.matrix(beta_0)))
            ##      beta_0_base <- Matrix::Matrix(vax_trans_red %x% beta_0)
            ##      beta_0_rtensor <- Matrix::Matrix(rTensor::kronecker_list(list(vax_trans_red, beta_0)))
        }

        beta_0 <- beta_0_fastmatrix ## beta_0_old

        ##   stopifnot(identical(beta_0_fastmatrix, beta_0_old))
        ##    stopifnot(identical(beta_0_base, beta_0_old))
        ##    stopifnot(identical(beta_0_rtensor, beta_0_old))

        ## expand names for beta_0 do colnames the same way in either case
        ## (whether beta_0 is a vector or matrix)
        col_names <- expand_names(original_col_names, vax_cat)
        ## prepare row names for vax-expanded beta_0 and update both rownames and
        ## colnames, depending on if beta_0
        ## vaxcats
        if (!is.null(original_row_names)) {
            row_names <- expand_names(original_row_names, vax_cat)
        } else {
            ## just update colnames (names) for vector
            row_names <- vax_cat
        }
        ## update dimnames for output
        dimnames(beta_0) <- list(row_names, col_names)
    } else {
        ## without vax
        if (do_variant(params)) {
            ## if we're doing variant with vax, we're just adjusting based on changes in transmissibility
            beta0 <- beta0 * calc_variant_adjustment(params, "unvax")
        }
    }

    ## assume that any matching values will be of the form "^%s_" where %s is something in Icats
    ## lapply(Icats, function(x) grep(sprintf("^%s_"), names(state))
    ## FIXME: we should be doing this by name, not assuming that all infectious compartments are expanded
    ##  into exactly 4 subcompartments, in order (but this should work for now??)
    if (has_testing(state = state)) { ## testified!
        if (has_age(params)) stop("can't combine age and testing yet")
        if (has_vax(params)) stop("can't combine vax and testing yet")
        beta_0 <- rep(beta_0, each = length(testcats))
        names(beta_0) <- unlist(lapply(Icats, function(x) paste0(x, testcats)))
        ## FIXME: also adjust _n, _p components?
        pos_vals <- grep("_t$", names(beta_0))
        beta_0[pos_vals] <- beta_0[pos_vals] * (1 - params[["iso_t"]])
    }
    if (!full) {
        return(beta_0)
    }

    ## By default, make a vector of zeroes for all the states,
    ## then fill in infectious ones
    ## without age and vax, just return vector
    if (!has_age(params) & !has_vax(params)) {
        beta <- setNames(numeric(length(state)), names(state))
        beta[names(beta_0)] <- beta_0
    } else {
        ## with age and vax, return a matrix
        beta <- matrix(0,
            nrow = nrow(beta_0), ncol = length(state),
            dimnames = list(rownames(beta_0), names(state))
        )
        beta[rownames(beta_0), colnames(beta_0)] <- matrix(beta_0)
    }
    return(beta)
}

## make_ratemat()
##' Create transition matrix
##'
##' Defines rates (per day) of flow \emph{from} compartment \code{i}
##' (row) \emph{to} compartment \code{j} (column).
##'
##' @details
##' The rates are as follows:
##'
##' \eqn{ S to E:  - (\beta_0 / N) S (C_a I_a + C_p I_p + (1-iso_m)C_m I_m + (1-iso_s)C_s I_s) }
##' \eqn{ E to I_a: }
##' \eqn{ E to I_p: }
##' \eqn{ ... }
##'
##' See \code{\link{read_params}} for parameter definitions.
##'
##' @note
##' Base version matches structure of Stanford/Georgia models
##' \itemize{
##'   \item flow diagram: see \url{http://covid-measures.stanford.edu/} 'model details' tab
##'         or \code{../pix/model_schematic.png}
##'   \item parameter definitions: see \code{params_CI_base.csv}, \code{params_ICU_diffs.csv}
##' }
##'
##' @param state named vector of states
##' @param params named vector of parameters
##' @param do_ICU include additional health utilization compartments
##' @param sparse return sparse matrix?
##' @param symbols return character (symbol) form? (FIXME: call adjust_symbols here rather than in show_ratemat()?)
##' @param indices return indices for lower-level stuff?
##' @return matrix of (daily) transition rates
##  *need* Matrix version of rowSums imported to handle sparse stuff below!!
##' @importFrom Matrix Matrix rowSums colSums
##' @examples
##' params <- read_params("ICU1.csv")
##' state <- make_state(params[["N"]],E0=params[["E0"]], use_eigvec=FALSE)
##' M <- make_ratemat(state,params)
##' if (require(Matrix)) {
##'    image(Matrix(M))
##' }
##' make_ratemat(state,params,symbols=TRUE)
##' @export
make_ratemat <- function(state, params, do_ICU = TRUE, sparse = FALSE,
                         symbols = FALSE, indices = FALSE) {
    pfun_opt <- getOption("macpan_pfun_method")
    if(has_vax(params) && (is.null(pfun_opt) || pfun_opt != "grep")) stop('please set `options(macpan_pfun_method = "grep")` at the top of your script')

    ## circumvent test code analyzers ... problematic ...
    alpha <- sigma <- gamma_a <- gamma_m <- gamma_s <- gamma_p <- NULL
    rho <- delta <- mu <- N <- E0 <- iso_m <- iso_s <- phi1 <- NULL
    phi2 <- psi1 <- psi2 <- psi3 <- c_prop <- c_delaymean <- c_delayCV <- NULL
    vax_response_rate <- vax_response_rate_R <- vax_alpha_dose1 <- NULL
    vax_mu_dose1 <- vax_alpha_dose2 <- vax_mu_dose2 <- NULL
    ## default values, will be masked (on purpose) by unpacking params/state
    nonhosp_mort <- 0
    ##
    np <- length(params)
    if (is.list(params)) {
        ## check param lengths (exclude setting-specific mistry params that will
        ## be of length 4,the number of settings, and not a multiple of age
        ## groups)
        nps <- lengths(params[!grepl("mistry_contact_rate_setting|mistry_fmats", names(params))])
        if (has_age(params)) {
            na <- length(attr(params, "age_cat"))
            bad_len <- which(!nps %in% c(1, na, na^2))
            if (length(bad_len) > 0) {
                stop(sprintf(
                    "elements of params must be length 1, %d or %d: %s",
                    na, na^2,
                    paste(names(params)[bad_len], collapse = ", ")
                ))
            }
        } else {
            if (!all(nps == 1)) stop("parameters should each be of length 1")
        }
    }
    state_names <- untestify_statenames(names(state))
    ns <- length(state_names)
    ## make param names locally available (similar to with())
    ## DON'T unpack states, we don't need them
    ## (the only state-dependent per capita rates are testing
    ## and infection, those get handled elsewhere)
    P <- as.list(params)
    unpack(P)
    ## blank matrix
    M <- Matrix::Matrix(0,
        nrow = ns, ncol = ns,
        dimnames = list(from = state_names, to = state_names)
    )
    if (symbols) {
      ## can't use a sparse matrix in this case ...
      dn <- dimnames(M)
      M <- as.matrix(M)
      storage.mode(M) <- "character"
      dimnames(M) <- dn
    }
    ## generic assignment function, indexes by regexp rather than string
    afun <- function(from, to, val) {
        if (!symbols) {
            M[pfun(from, to, M)] <<- val
        } else {
            M[pfun(from, to, M)] <<- deparse(substitute(val))
        }
    }

    ## fill entries in ratemat

    ## calculate FOI
    beta_array <- make_beta(state, params)
    ## FIXME: call update_foi() here?
    if ("numeric" %in% class(beta_array)) {
        ## dot product of beta_array (vec) with states (I classes) to get FOI
        afun("S", "E", sum(beta_array * state[names(beta_array)]))
    } else {
        ## dot product of each row of beta_array (corresponding to a different S subcategory) with states (I classes) to get FOIs (plural, one per susceptible class)
        afun("S", "E", beta_array %*% state[colnames(beta_array)])
    }

    ## fill other parameters
    afun("E", "Ia", alpha * sigma)
    afun("E", "Ip", (1 - alpha) * sigma)
    afun("Ia", "R", gamma_a)
    afun("Ip", "Im", mu * gamma_p)
    afun("Ip", "Is", (1 - mu) * gamma_p)
    afun("Im", "R", gamma_m)

    ## fill hospital-related parameters, depending on ICU model selected
    if (!do_ICU) {
        ## simple hospital model as in Stanford/CEID
        afun("Is", "H", gamma_s)
        afun("H", "D", delta * rho)
        afun("H", "R", (1 - delta) * rho)
    } else {
        ## FIXME: A better term than "acute" to mean the opposite of intensive?
        ## four-way split (direct to D, acute care, ICU/survive, ICUD/die)?
        afun("Is", "H", (1 - nonhosp_mort) * phi1 * gamma_s)
        afun("Is", "ICUs", (1 - nonhosp_mort) * (1 - phi1) * (1 - phi2) * gamma_s)
        afun("Is", "ICUd", (1 - nonhosp_mort) * (1 - phi1) * phi2 * gamma_s)
        afun("Is", "D", nonhosp_mort * gamma_s)
        afun("ICUs", "H2", psi1) ## ICU to post-ICU acute care
        afun("ICUd", "D", psi2) ## ICU to death
        afun("H2", "R", psi3) ## post-ICU to discharge
        ## H now means 'acute care' only; all H survive & are discharged
        afun("H", "D", 0)
        afun("H", "R", rho) ## all acute-care survive
        if (any(grepl("^X", colnames(M)))) {
            ## FIXME: check for age?
            afun("Is", "X", M[pfun("Is", "H", M)]) ## assuming that hosp admissions mean *all* (acute-care + ICU)
        }
    }

    ## FIXME: adjust epi parameters for vaxprotect strata here (e.g. increase asymptomatic proportion, alpha, and set proportion of symptomatic infections that are mild, mu, to 1, as well as non-hosp mortality to 0)

    ## vaccination-related rates
    if (has_vax(params)) {
        ## get vax categories
        vax_cat <- get_vax(params)
        model_type <- get_vax_model_type(vax_cat)

        ## add vaccine allocation rates (from unvax to vaxdose)
        M <- add_updated_vaxrate(state, params, M)

        ## add vaccine immune response rate for not active infections
        ## dose1 -> protect1
        afun(
            paste0("S_.*", vax_cat[2]),
            paste0("S_.*", vax_cat[3]),
            vax_response_rate
        )

        afun(
            paste0("R_.*", vax_cat[2]),
            paste0("R_.*", vax_cat[3]),
            vax_response_rate_R
        )

        ## dose2 -> protect2
        if (model_type == "twodose") {
            afun(
                paste0("S_.*", vax_cat[4]),
                paste0("S_.*", vax_cat[5]),
                vax_response_rate
            )

            afun(
                paste0("R_.*", vax_cat[4]),
                paste0("R_.*", vax_cat[5]),
                vax_response_rate_R
            )
        }

        ## modify epidemiological parameters for vaxprotect1
        ## individuals
        afun(
            paste0("E_.*", vax_cat[3]),
            paste0("Ia_.*", vax_cat[3]), vax_alpha_dose1 * sigma
        )
        afun(
            paste0("E_.*", vax_cat[3]),
            paste0("Ip_.*", vax_cat[3]),
            (1 - vax_alpha_dose1) * sigma
        )
        afun(
            paste0("Ip_.*", vax_cat[3]),
            paste0("Im_.*", vax_cat[3]),
            vax_mu_dose1 * gamma_p
        )
        afun(
            paste0("Ip_.*", vax_cat[3]),
            paste0("Is_.*", vax_cat[3]),
            (1 - vax_mu_dose1) * gamma_p
        )

        if (model_type == "twodose") {
            ## carry over epi param adjustments from vaxprotect1
            ## to vaxdose2 layer
            afun(
                paste0("E_.*", vax_cat[4]),
                paste0("Ia_.*", vax_cat[4]), vax_alpha_dose1 * sigma
            )
            afun(
                paste0("E_.*", vax_cat[4]),
                paste0("Ip_.*", vax_cat[4]),
                (1 - vax_alpha_dose1) * sigma
            )
            afun(
                paste0("Ip_.*", vax_cat[4]),
                paste0("Im_.*", vax_cat[4]),
                vax_mu_dose1 * gamma_p
            )
            afun(
                paste0("Ip_.*", vax_cat[4]),
                paste0("Is_.*", vax_cat[4]),
                (1 - vax_mu_dose1) * gamma_p
            )

            ## adjust epi parameters for fully-vaccinated people
            afun(
                paste0("E_.*", vax_cat[5]),
                paste0("Ia_.*", vax_cat[5]), vax_alpha_dose2 * sigma
            )
            afun(
                paste0("E_.*", vax_cat[5]),
                paste0("Ip_.*", vax_cat[5]),
                (1 - vax_alpha_dose2) * sigma
            )
            afun(
                paste0("Ip_.*", vax_cat[5]),
                paste0("Im_.*", vax_cat[5]),
                vax_mu_dose2 * gamma_p
            )
            afun(
                paste0("Ip_.*", vax_cat[5]),
                paste0("Is_.*", vax_cat[5]),
                (1 - vax_mu_dose2) * gamma_p
            )
        }
    }

    if (sparse) {
        M <- Matrix::Matrix(M)
    } else {
        M <- as.matrix(M)
    }

    return(M)
}

##' calculate only updated force of infection
##' at present, this is the only state-dependent \emph{per capita} rate
##' maybe more efficient than modifying & returning the whole matrix
##' @inheritParams make_ratemat
##' @param beta vector or matrix of transmission rates, where (length(beta) | ncol(beta)) == length(state)
##' @export
## FIXME DRY from make_ratemat
update_foi <- function(state, params, beta) {

    ## update infection rate
    if (is.matrix(beta)) {
        ## e.g. beta is a matrix when the transmission parameters are age-structured
        if (length(state) != ncol(beta)) stop("number of columns of beta must match length of state vector")
        foi <- beta %*% state
    } else {
        if (length(state) != length(beta)) {
            stop("length of state and beta are not the same")
        }
        foi <- sum(state * beta)
    }
    if (has_zeta(params)) {
        if (has_age(params)) stop("phenomenological heterogeneity (zeta != 0) untested with age-structed params")
        ## suppose zeta is a vector zeta1, zeta2, zeta3, ...
        ##  we also need parameters   zeta_break1, zeta_break2 (0<zbi<1)
        ##  one *fewer* break parameter than zeta_i value
        ## if 0< S/N < zeta_break1   -> zeta1
        ##  zeta_break1 < S/N < zeta_break2 -> zeta2
        ## ...
        ##  zeta_breakx < S/N < 1  -> zetax
        Susc_frac <- 1 / params[["N"]] * sum(state[grep("^S_?", names(state))])
        if (any(grepl("zeta[0-9]", names(params)))) {
            zeta <- with(
                as.list(params),
                if (Susc_frac < zeta_break) zeta1 else zeta2
            )
        }
        ## alternately could just make it a vector
        ## ... but this messes with age-structured stuff
        ## if (length(zeta)>0) {
        ## ...
        ## }
        ## ???? het at pop level or sub-category level??
        foi <- foi * with(as.list(params), Susc_frac^zeta)
    }
    return(foi)
}

## update the entire rate matrix; we need this when we are doing testify models
##  because we expect testing rates to change daily ...
update_ratemat <- function(ratemat, state, params, testwt_scale = "N") {
    if (inherits(ratemat, "Matrix") && has_testing(params)) {
        aa <- c("wtsvec", "posvec", "testing_time")
        saved_attrs <- setNames(lapply(aa, attr, x = ratemat), aa)
    }
    ## update testing flows. DO THIS FIRST, before updating foi: **** assignment via cbind() to Matrix objects loses attributes???
    if (has_testing(state)) {
        testing_time <- attr(ratemat, "testing_time") ## ugh  (see **** above)
        ## positions of untested, positive-waiting, negative-waiting compartments
        ## (flows from _n, _p to _t, or back to _u, are always at per capita rate omega, don't need
        ##  to be updated for changes in state)
        ## FIXME: backport to testify?
        u_pos <- grep("_u$", rownames(ratemat))
        p_pos <- grep("_p$", rownames(ratemat))
        n_pos <- grep("_n$", rownames(ratemat))
        ## original/unscaled prob of positive test by compartment, testing weights by compartment
        posvec <- attr(ratemat, "posvec")
        if (is.null(posvec)) stop("expected ratemat to have a posvec attribute")
        wtsvec <- attr(ratemat, "wtsvec")
        if (is.null(wtsvec)) stop("expected ratemat to have a wtsvec attribute")
        ## scaling ...
        testing_intensity <- params[["testing_intensity"]]
        testing_tau <- params[["testing_tau"]]
        N0 <- params[["N"]]
        W <- sum(wtsvec * state[u_pos])
        sc <- switch(testwt_scale,
            none = 1,
            N = N0 / W,
            sum_u = sum(state[u_pos]) / W,
            sum_smooth = {
                rho <- testing_intensity
                tau <- testing_tau
                tau * N0 / (tau * W + rho * N0)
                ## NOTE 'smoothing' doc has numerator rho*tau*N0,
                ## but testing intensity (rho) is included in ratemat
                ## calculation below ...
            }
        )
        ratemat[cbind(u_pos, n_pos)] <- testing_intensity * sc * wtsvec * (1 - posvec)
        ratemat[cbind(u_pos, p_pos)] <- testing_intensity * sc * wtsvec * posvec
        if (testing_time == "sample") {
            N_pos <- which(rownames(ratemat) == "N")
            P_pos <- which(rownames(ratemat) == "P")
            ratemat[cbind(u_pos, N_pos)] <- ratemat[cbind(u_pos, n_pos)]
            ratemat[cbind(u_pos, P_pos)] <- ratemat[cbind(u_pos, p_pos)]
        }
    }

    ## update FOIs
    ratemat[pfun("S", "E", ratemat)] <- update_foi(state, params, make_beta(state, params))

    ## update vaccine allocation rates
    if (has_vax(params)) {
        ratemat <- add_updated_vaxrate(state, params, ratemat)
    }

    ## ugh, restore attributes if necessary
    if (inherits(ratemat, "Matrix") && has_testing(params)) {
        for (a in aa) {
            attr(ratemat, a) <- saved_attrs[[a]]
        }
    }
    return(ratemat)
}

## do_step()
##' Take a single simulation time step
##' @inheritParams make_ratemat
##' @param ratemat transition matrix
##' @param dt time step (days)
##' @param do_hazard use hazard calculation?
##' @param do_exponential prevent outflow of susceptibles, to create a pure-exponential process?
##' @param stoch_proc stochastic process error?
##' @param testwt_scale how to scale testing weights? "none"=use original weights as specified;
##' "N" = multiply by (pop size)/(sum(wts*state[u_pop])); "sum_u" = multiply by (sum(state[u_pop])/(sum(wts*state[u_pop])))
##' @export
##' @examples
##' params1 <- read_params("ICU1.csv")
##' state1 <- make_state(params=params1)
##' M <- make_ratemat(params=params1, state=state1)
##' s1A <- do_step(state1,params1, M, stoch_proc=TRUE)
do_step <- function(state, params, ratemat, dt = 1,
                    do_hazard = TRUE,
                    stoch_proc = FALSE,
                    do_exponential = FALSE,
                    testwt_scale = "N") {
    x_states <- c("X", "N", "P") ## weird parallel accumulators
    ## if vax accumulator is in the state vec, add it to the list of parallel accumulators
    if ("V" %in% attr(state, "epi_cat")) {
        x_states <- c(x_states, "V")
    }
    p_states <- exclude_states(names(state), x_states)
    ## FIXME: check (here or elsewhere) for non-integer state and process stoch?
    ## cat("do_step beta0",params[["beta0"]],"\n")
    ratemat <- update_ratemat(ratemat, state, params, testwt_scale = testwt_scale)
    if (!stoch_proc || (!is.null(s <- params[["proc_disp"]]) && s < 0)) {
        if (!do_hazard) {
            ## from per capita rates to absolute changes
            flows <- col_multiply(ratemat, state) * dt ## sweep(ratemat, state, MARGIN=1, FUN="*")*dt
        } else {
            ## FIXME: change var names? {S,E} is a little confusing (sum, exp not susc/exposed)
            ## use hazard function: assumes exponential change
            ## (constant per capita flows) rather than linear change
            ## (constant absolute flows) within time steps
            ## handle this as in pomp::reulermultinom,
            ## i.e.
            ##    S = sum(r_i)   ## total rate
            ##    p_{ij}=(1-exp(-S*dt))*r_j/S
            ##    p_{ii}= exp(-S*dt)
            S <- rowSums(ratemat)
            E <- exp(-S * dt)
            ## prevent division-by-0 (boxes with no outflow) problems (FIXME: DOUBLE-CHECK)
            norm_sum <- ifelse(S == 0, 0, state / S)
            flows <- (1 - E) * col_multiply(ratemat, norm_sum) ## sweep(ratemat, norm_sum, MARGIN=1, FUN="*")
            diag(flows) <- 0 ## no flow
        }
    } else {
        flows <- ratemat ## copy structure
        flows[] <- 0
        for (i in seq(length(state))) {
            ## FIXME: allow Dirichlet-multinomial ?
            dW <- dt
            if (!is.na(proc_disp <- params[["proc_disp"]])) {
                dW <- pomp::rgammawn(sigma = proc_disp, dt = dt)
            }
            ## FIXME: need to adjust for non-conserving accumulators!
            flows[i, -i] <- pomp::reulermultinom(
                n = 1,
                size = state[[i]],
                rate = ratemat[i, -i],
                dt = dW
            )
        }
    }

    if (!do_exponential) {
        outflow <- rowSums(flows[, p_states])
    } else {
        ## want to zero out outflows from S to non-S compartments
        ##  (but leave the inflows - thus we can't just zero these flows
        ##   out in the rate matrix!)
        S_pos <- grep("^S", rownames(ratemat), value = TRUE)
        notS_pos <- grep("^[^S]", colnames(ratemat), value = TRUE)
        notS_pos <- setdiff(notS_pos, x_states)
        outflow <- setNames(numeric(ncol(flows)), colnames(flows))
        ## only count outflows from S_pos to other S_pos (e.g. testing flows)
        outflow[S_pos] <- rowSums(flows[S_pos, S_pos, drop = FALSE])
        ## count flows from infected etc. to p_states (i.e. states that are *not* parallel accumulators)
        outflow[notS_pos] <- rowSums(flows[notS_pos, p_states])
    }
    inflow <- colSums(flows)
    state <- state - outflow + inflow
    ## check conservation (*don't* check if we are doing an exponential sim, where we
    ##  allow infecteds to increase without depleting S ...)
    MP_badsum_action <- getOption("MP_badsum_action", "warning")
    MP_badsum_tol <- getOption("MP_badsum_tol", 1e-12)
    if (!do_exponential &&
        !(MP_badsum_action == "ignore") &&
        !stoch_proc) ## temporary: adjust reulermultinom to allow for x_states ...
        {
            calc_N <- sum(state[p_states])
            if (!isTRUE(all.equal(calc_N, sum(params[["N"]]), tolerance = MP_badsum_tol))) {
                msg <- sprintf("sum(states) != original N (delta=%1.2g)", sum(params[["N"]]) - calc_N)
                get(MP_badsum_action)(msg)
            }
        } ## not exponential run or stoch proc or ignore-sum
    return(state)
    ## Why is this throwing warnings in make_state?
    ##    if(any(state < -sqrt(.Machine$double.eps))){
    ##        warning('End of run_sim_range check: One or more state variables is negative, below -sqrt(.Machine$double.eps)')
    ##    }
    ##    else if(any(state < 0)){
    ##        warning('End of run_sim_range check: One or more state variables is negative, below -sqrt(.Machine$double.eps) and 0.')
    ##    }
}

## global flag
deprecate_timepars_warning <- FALSE

## run_sim()
##' Run pandemic simulation
##' @inheritParams do_step
##' @param start_date starting date (Date or character, any sensible D-M-Y format)
##' @param end_date ending date (ditto)
##' @param params_timevar four-column data frame containing columns 'Date'; 'Symbol' (parameter name/symbol); 'Value'; and 'Type' (\code{rel_orig} (relative to value at time zero); \code{rel_prev} (relative to previous value); or \code{abs} (absolute value at specified time). If a 'Relative_value' column is present, it is renamed to 'Value' and 'Type' is set to \code{rel_orig} for back-compatibility.
##' @param dt time step for \code{\link{do_step}}
##' @param ratemat_args additional arguments to pass to \code{\link{make_ratemat}}
##' @param step_args additional arguments to pass to \code{\link{do_step}}
##' @param ndt number of internal time steps per time step
##' @param stoch a logical vector with elements "obs" (add obs error?) and "proc" (add process noise?)
##' @param stoch_start dates on which to enable stochasticity (vector of dates with names 'proc' and 'obs')
##' @param condense if \code{TRUE}, use \code{\link{condense.pansim}} to reduce the number of variables in the output (in particular, collapse subclasses and return only one \code{I}, \code{H}, and \code{ICU} variable)
##' @param condense_args arguments to pass to \code{\link{condense}} (before adding observation error)
##' @param use_ode integrate via ODE rather than discrete step?
##' @param ode_args additional arguments to \code{\link[deSolve]{ode}}
##' @param use_flex use \code{flexmodel} approach (experimental)
##' @examples
##' params <- read_params("ICU1.csv")
##' paramsS <- update(params,c(proc_disp=0.1,obs_disp=100))
##' paramsSz <- update(paramsS, zeta=5)
##' state <- make_state(params=params)
##' time_pars <- data.frame(Date=c("2020-03-20","2020-03-25"),
##'                        Symbol=c("beta0","beta0"),
##'                        Relative_value=c(0.7,0.1),
##'                        stringsAsFactors=FALSE)
##' res1 <- run_sim(params,state,start_date="2020-02-01",end_date="2020-06-01")
##' res1X <- run_sim(params,state,start_date="2020-02-01",end_date="2020-06-01",
##'                  condense_args=list(keep_all=TRUE))
##' res1_S <- update(res1, params=paramsS, stoch=c(obs=TRUE, proc=TRUE))
##' res1_t <- update(res1, params_timevar=time_pars)
##' res1_S_t <- update(res1_S, params_timevar=time_pars)
##' res2_S_t <- update(res1_S_t,params=update(paramsS, proc_disp=0.5))
##' res3_S_t <- update(res2_S_t,stoch_start="2020-04-01")
##' res3_Sz <- update(res1_S, params=paramsSz)
##' plot(res3_Sz,log=TRUE,log_lwr=1e-4)
##' @importFrom stats rnbinom na.exclude napredict
##' @param verbose print messages (e.g. about time-varying parameters)?
##' @export
## FIXME: automate state construction better
run_sim <- function(params,
                    state = NULL,
                    start_date = "2020-03-20",
                    end_date = "2020-05-1",
                    params_timevar = NULL,
                    dt = 1,
                    ndt = 1 ## FIXME: change default after testing?
                    , stoch = c(obs = FALSE, proc = FALSE),
                    stoch_start = NULL,
                    ratemat_args = NULL,
                    step_args = list(),
                    ode_args = list(),
                    use_ode = FALSE,
                    condense = TRUE,
                    condense_args = NULL,
                    verbose = FALSE,
                    use_flex = FALSE,
                    flexmodel = NULL,
                    obj_fun = NULL) {

  if (!is.null(flexmodel) | !is.null(obj_fun)) use_flex = TRUE
  if (use_flex) {
    ## tmb/c++ computational approach (experimental)
    ## (https://canmod.net/misc/flex_specs)

    start_date <- as.Date(start_date)
    end_date <- as.Date(end_date)
    if(!is.null(params_timevar)) params_timevar$Date = as.Date(params_timevar$Date)

    # update parameters (or the whole model structure)
    # -- important in calibration situations where the parameters
    #    are changing each iteration of the optimizer
    if (!is.null(flexmodel)) {
      flexmodel$params = expand_params_S0(params, 1-1e-5)
      if (spec_ver_gt('0.1.0') & (isTRUE(nrow(flexmodel$timevar$piece_wise$schedule) > 0))) {
        if (is.null(params_timevar)) {
          params_timevar = flexmodel$timevar$piece_wise$schedule
        }
        # FIXME: inefficient brute-force reordering
        s = arrange(flexmodel$timevar$piece_wise$schedule, Symbol, Date)
        ptv = arrange(params_timevar, Symbol, Date)
        s$last_tv_mult = ptv$Value
        flexmodel$timevar$piece_wise$schedule = arrange(s, breaks, tv_spi)
      }
    } else {
      # FIXME: can't assume base model
      #        -- at least throw error if the model has structure
      flexmodel <- make_base_model(
        params = params,
        state = state,
        start_date = start_date,
        end_date = end_date,
        params_timevar = params_timevar,
        step_args = step_args
      )
    }

    if (is.null(state)) {
      state = flexmodel$state
    }

    if(!(spec_ver_gt('0.1.0') & isTRUE(flexmodel$do_make_state)) | is.null(obj_fun)) {
      # if not making the state on the c++ side, need to make
      # a new tmb fun with the new state that has come in
      # from up the call stack
      flexmodel$state = state
      obj_fun <- tmb_fun(flexmodel)
    }

    ## simulate trajectories based on new parameters
    full_param_vec = tmb_params(flexmodel)
    tmb_sims <- obj_fun$simulate(full_param_vec)
    res = (tmb_sims
      %>% getElement("concatenated_state_vector")
      %>% structure_state_vector(flexmodel$iters, names(flexmodel$state))
    )
    res <- dfs(date = seq(flexmodel$start_date, flexmodel$end_date, by = 1), res)

    # look for foi -- TODO: formalize the definition of foi so that we
    #                       can reliably check -- this is too model-specific
    foi = (flexmodel$tmb_indices$updateidx
      %>% names
      %>% strsplit("_to_")
      %>% lapply(function(x) {
        # vax_cat may not exist and therefore be a blank string
        vax_cat = sub("^E(_|$)", "", x[2])
        setNames(
          startsWith(x[1], "S") & startsWith(x[2], "E"),
          ifelse(vax_cat == '', 'foi', 'foi' %_% vax_cat))
      })
      %>% unlist
      %>% which
      %>% lapply(function(foi_off) {
        tmb_sims$concatenated_ratemat_nonzeros[
          seq(from = foi_off,
              by = length(flexmodel$tmb_indices$updateidx),
              length.out = flexmodel$iters + 1)
        ]
      })
      %>% as.data.frame
    )

    params0 <- params
    state0 = state
    if (spec_ver_gt('0.1.0') & isTRUE(flexmodel$do_make_state)) {
      # if the initial state came from tmb ...
      state0[] = tmb_sims$concatenated_state_vector[seq_along(state)]
    }

    res = cbind(res, foi)
  } else {

    call <- match.call()

    if (is.na(sum(params[["N"]]))) stop("no population size specified; set params[['N']]")

    ## FIXME: *_args approach (specifying arguments to pass through to
    ##  make_ratemat() and do_step) avoids cluttering the argument
    ##  list, but may be harder to translate to lower-level code
    if (dt != 1) warning("nothing has been tested with dt!=1")
    start_date <- as.Date(start_date)
    end_date <- as.Date(end_date)
    if (!is.null(stoch_start)) {
      stoch_start <- setNames(as.Date(stoch_start), names(stoch_start))
    }
    if (length(stoch_start) == 1) stoch_start <- c(obs = stoch_start, proc = stoch_start)
    date_vec <- seq(start_date, end_date, by = dt)
    nt <- length(date_vec)
    step_args <- c(step_args, list(stoch_proc = stoch[["proc"]]))
    drop_last <- function(x) {
      x[seq(nrow(x) - 1), ]
    }
    if (is.null(state)) state <- make_state(params = params, testify = FALSE)
    if (has_testing(params = params)) {
      state <- expand_stateval_testing(state, params = params)
    }

    if (has_vax(params) && ndt != 1) stop("the model with vaccination only works with ndt = 1")

    ## thin wrapper: we may have to recompute this later and don't want to
    ##  repeat both make_ratemat() and all of the testify stuff ...
    ## (1) is all of this idempotent, i.e. will re-running expand_stateval_testing break anything if unnecessary?
    ## (2) does this imply that testify/testing stuff should be refactored/go elsewhere?
    ## (3) this is reminiscent of the expand/don't-expand hoops that we go through in the eigenvector/state calculation
    make_M <- function() {
      M <- make_ratemat(state = state, params = params)
      if (has_testing(params = params)) {
        if (!is.null(ratemat_args$testify)) {
          warning("'testify' no longer needs to be passed in ratemat_args")
        }
        testing_time <- ratemat_args$testing_time
        if (is.null(testing_time)) {
          warning("setting testing time to 'sample'")
          testing_time <- "sample"
        }
        M <- testify(M, params, testing_time = testing_time)
      }
      return(M)
    }
    M <- make_M()

    state0 <- state
    params0 <- params ## save baseline (time-0) values
    ## no explicit switches, and (no process error) or (process error for full time);
    ## we will be able to run the whole sim directly
    if (is.null(params_timevar) && (!stoch[["proc"]] || is.null(stoch_start))) {
      switch_times <- NULL
    } else {
      if (is.null(params_timevar)) {
        ## starting times for process/obs error specified, but no other time-varying parameters;
        ##  we need an empty data frame with the right structure so we can append the process-error switch times
        params_timevar <- dfs(Date = as.Date(character(0)), Symbol = character(0), Value = numeric(0), Type = character(0))
      } else {
        ## check column names
        names(params_timevar) <- capitalize(names(params_timevar))
        if (identical(
          names(params_timevar),
          c("Date", "Symbol", "Relative_value")
        )) {
          if (!deprecate_timepars_warning) {
            warning("specifying params_timevar with Relative_value is deprecated: auto-converting (reported once per session)")
            deprecate_timepars_warning <- TRUE
          }
          names(params_timevar)[3] <- "Value"
          params_timevar <- data.frame(params_timevar, Type = "rel_orig")
        }
        npt <- names(params_timevar)
        if (!identical(
          npt,
          c("Date", "Symbol", "Value", "Type")
        )) {
          stop("params_timevar: has wrong names: ", paste(npt, collapse = ", "))
        }
        params_timevar$Date <- as.Date(params_timevar$Date)
        ## tryCatch(
        ##     params_timevar$Date <- as.Date(params_timevar$Date),
        ##     error=function(e) stop("Date column of params_timevar must be a Date, or convertible via as.Date"))
        params_timevar <- params_timevar[order(params_timevar$Date), ]
      }
      ## append process-observation switch to timevar
      if (stoch[["proc"]] && !is.null(stoch_start)) {
        params_timevar <- rbind(
          params_timevar,
          dfs(
            Date = stoch_start[["proc"]],
            Symbol = "proc_disp", Value = 1, Type = "rel_orig"
          )
        )
        params[["proc_disp"]] <- -1 ## special value: signal no proc error
      }
      switch_dates <- params_timevar[["Date"]]
      ## match specified times with time sequence
      switch_times <- match(switch_dates, date_vec)
      if (any(is.na(switch_times))) {
        bad <- which(is.na(switch_times))
        stop("non-matching dates in params_timevar: ", paste(switch_dates[bad], collapse = ","))
      }
      if (any(switch_times == length(date_vec))) {
        ## drop switch times on final day
        warning("dropped switch times on final day")
        switch_times <- switch_times[switch_times < length(date_vec)]
      }
    } ## steps

    if (is.null(switch_times)) {
      res <- thin(
        ndt = ndt,
        do.call(
          run_sim_range,
          nlist(params,
                state,
                nt = nt * ndt,
                dt = dt / ndt,
                M,
                use_ode,
                ratemat_args,
                step_args
          )
        )
      )
    } else {
      t_cur <- 1
        ## want to *include* end date
        switch_times <- switch_times + 1
        ## add beginning and ending time
        times <- c(1, unique(switch_times), nt + 1)
        resList <- list()
        ## for switch time indices
        for (i in seq(length(times) - 1)) {
            recompute_M <- FALSE
            for (j in which(switch_times == times[i])) {
                ## reset all changing params
                s <- params_timevar[j, "Symbol"]
                v <- params_timevar[j, "Value"]
                t <- params_timevar[j, "Type"]
                if (t == "abs") {
                  if (length(unique(params[[s]])) > 1) {
                    stop("attempting to replace a vector-valued parameter with a scalar value")
                  }
                  params[[s]] <- v
                } else {
                  old_param <- switch(t,
                      ## this should work even if params0[[s]] is a vector
                      rel_orig = params0[[s]],
                      rel_prev = params[[s]],
                      stop("unknown time_params type ", t)
                  )
                  params[[s]] <- old_param * v
                }
                if (s == "proc_disp") {
                    state <- round(state)
                }
                if (verbose) {
                    cat(sprintf(
                        "changing value of %s from original %f to %f at time step %d\n",
                        s, params0[[s]], params[[s]], i
                    ))
                }

                  if (!s %in% "beta0") { ## also testing rates?
                      recompute_M <- TRUE
                  }
                  ## FIXME: so far still assuming that params only change foi
                  ## if we change another parameter we will have to recompute M
              }
              if (recompute_M) {
                  M <- make_M()
              }

              resList[[i]] <- drop_last(
                  thin(
                      ndt = ndt,
                      do.call(
                          run_sim_range,
                          nlist(params,
                              state,
                              nt = (times[i + 1] - times[i] + 1) * ndt,
                              dt = dt / ndt,
                              M,
                              use_ode,
                              ratemat_args,
                              step_args,
                              ode_args
                          )
                      )
                  )
              )
              state <- attr(resList[[i]], "state")
              t_cur <- times[i]
          }
          ## combine
          res <- do.call(rbind, resList)
          ## add last row
          ## res <- rbind(res, attr(resList[[length(resList)]],"state"))
      }
      ## drop internal stuff
      ## res <- res[,setdiff(names(res),c("t","foi"))]
      res <- dfs(date = seq(start_date, end_date, by = dt), res)
      res <- res[, !names(res) %in% "t"] ## we never want the internal time vector
    } # not use_flex

    ## condense here
    if (condense) {
      ## FIXME: will it break anything if we call condense with 'params'
      ## (modified by time-varying stuff) instead of params0? Let's find out!
      res <- do.call(condense.pansim, c(
          list(res,
              params = params0,
              cum_reports = FALSE,
              het_S = has_zeta(params0)
          ),
          condense_args
      ))
    }
    if (stoch[["obs"]]) {
        if (has_zeta(params)) params[["obs_disp_hetS"]] <- NA ## hard-code skipping obs noise
        ## do observation error here
        ## FIXME: warn on mu<0 ? (annoying with ESS machinery)
        m <- res[, -1] ## drop time component
        if (!is.null(stoch_start)) {
            ## only add stochastic obs error after specified date
            m <- m[res$date > stoch_start[["obs"]], ]
        }
        m_rows <- nrow(m)
        for (i in seq(ncol(m))) {
            nm <- names(m)[i]
            if ((vn <- paste0("obs_disp_", nm)) %in% names(params)) {
                ## variable-specific dispersion specified
                d <- params[[vn]]
            } else {
                d <- params[["obs_disp"]]
            }
            if (!is.na(d)) {
                ## rnbinom, skipping NA values (as below)
                m[[i]] <- suppressWarnings(rnbinom(m_rows, mu = m[[i]], size = d))
            }
        }
        res[seq(nrow(res) - m_rows + 1, nrow(res)), -1] <- m
    }
    ## add cum reports *after* adding obs error
    ## TODO: how do reports and deaths interact with flexmodel/tmb approach?
    if ("report" %in% names(res)) res$cumRep <- cumsum(ifelse(!is.na(unname(unlist(res$report))), unname(unlist(res$report)), 0))
    if ("death" %in% names(res)) res$D <- cumsum(ifelse(!is.na(res$death), res$death, 0))

    ## repair age groups in result names (if present)
    res <- repair_names_age(res)

    ## store everything as attributes
    attr(res, "params") <- params0
    attr(res, "state0") <- state0
    attr(res, "start_date") <- start_date
    attr(res, "end_date") <- end_date
    attr(res, "call") <- call
    if(use_flex) {
      attr(res, "flexmodel") <- flexmodel
      attr(res, "ad_fun") <- obj_fun

    }
    attr(res, "params_timevar") <- params_timevar
    ## attr(res,"final_state") <- state
    class(res) <- c("pansim", "data.frame")

    state_names <- names(res)
    state_names_indices <- first_letter_cap(state_names)
    state_names <- state_names[state_names_indices]

    ## FIXME: why do we need to restrict to res[state_names] ?
    ## because e.g. report can be NaN for the first part of the sim (by default, the first 15 days, set in start_date_offset)
    if (isTRUE(any(res[state_names] < -sqrt(.Machine$double.eps)))) {
        state_vars <- (res %>% select(state_names))
        below_zero_lines <- (rowSums(state_vars < -sqrt(.Machine$double.eps)) > 0)

        warning(
            "End of run_sim check: One or more state variables is negative at some time, below -sqrt(.Machine$double.eps). Check following message for details \n ",
            paste(utils::capture.output(print(res[below_zero_lines, ])), collapse = "\n")
        )
    }
    else if (isTRUE(any(res[state_names] < 0))) {
        state_vars <- (res %>% select(state_names))
        below_zero_lines <- (rowSums(state_vars < 0) > 0)

        warning(
            "End of run_sim check: One or more state variables is negative at some time, between -sqrt(.Machine$double.eps) and 0. Check following message for details \n ",
            paste(utils::capture.output(print(res[below_zero_lines, ])), collapse = "\n")
        )
    }
    return(res)
}



##' retrieve parameters from a CSV file
##'
##' @details
##' The parameters that must be set are:
##'
##' \eqn{   N:  }  population size
##'
##' \eqn{   \beta_0:  }  transmission rate
##'
##' \eqn{   1/\sigma:  }  mean \emph{latent} period
##'
##' \eqn{   1/\gamma_a:  }  mean \emph{infectious} period for asymptomatic individuals
##'
##' \eqn{   ... }
##'

##' generate initial state vector
##' @param N population size
##' @param E0 initial number exposed
##' @param type (character) specify what model type this is intended
##'     for (e.g., \code{"ICU1"}, \code{"CI"}); determines state names
##' @param state_names vector of state names, must include S and E
##' @param params parameter vector (looked in for N and E0)
##' @param x proposed (named) state vector; missing values will be set
##' @param normalize (logical) should the state vector be normalized to sum to 1?
##' @param use_eigvec use dominant eigenvector to distribute non-Susc values
##' to zero: default is to set this to \code{TRUE} if \code{params} is non-NULL
##' @param testify expand state vector to include testing compartments (untested, neg waiting, pos waiting, pos received) ?
##' @param ageify expand state vector to include different age groups (??)
##' @param vaxify expand state vector to include groups that have 1-2 vaccine doses (??)
##' @note \code{"CI"} refers to the Stanford group's
##'     "covid intervention" model.
##' @export
##' @examples
##' p <- read_params("ICU1.csv")
##' make_state(N=1e6,E0=1)
##' make_state(params=p)
##  FIXME: can pass x, have a name check, fill in zero values
make_state <- function(N = params[["N"]],
                       E0 = params[["E0"]],
                       type = "ICU1h",
                       state_names = NULL,
                       use_eigvec = NULL,
                       params = NULL,
                       x = NULL,
                       normalize = FALSE,
                       ageify = NULL,
                       vaxify = NULL,
                       testify = NULL) {
    if (is.null(ageify)) ageify <- (!is.null(params) && has_age(params)) || (!is.null(x) && any(grepl("\\+$", names(x))))
    if (is.null(vaxify)) vaxify <- (!is.null(params) && has_vax(params)) || (!is.null(x) && any(grepl("vax", names(x))))
    if (is.null(testify)) testify <- !is.null(params) && has_testing(params = params)
    ## error if use_eigvec was **explicitly requested** (equiv !missing(use_eigvec)) && no params
    if (isTRUE(use_eigvec) && is.null(params)) stop("must specify params")
    if (is.null(use_eigvec)) use_eigvec <- !is.null(params)

    ## select vector of epi state names
    state_names <- switch(type,
        ## "X" is a hospital-accumulator compartment (diff(X) -> hosp)
        ## "V" is the vaccination accumulator compartment
        ICU1h = c("S", "E", "Ia", "Ip", "Im", "Is", "H", "H2", "ICUs", "ICUd", "D", "R", "X", "V"),
        ICU1 = c("S", "E", "Ia", "Ip", "Im", "Is", "H", "H2", "ICUs", "ICUd", "D", "R"),
        CI =   c("S", "E", "Ia", "Ip", "Im", "Is", "H", "D", "R"),
        ## Add a test case which should result in a thrown warning
        test_warning_throw = c("s", "e", "Ia", "Ip", "Im", "Is", "H", "H2", "ICUs", "ICUd", "d", "R", "X"),
        stop("unknown type")
    )
    epi_cat <- state_names

    ## set up output state vector, depending on what strata are requested
    state <- setNames(numeric(length(state_names)), state_names)
    if (ageify) {
        if (!is.null(params)) state <- expand_state_age(state, attr(params, "age_cat"))
        if (!is.null(x)) {
            state_names <- names(x) ## update state names based on those provided in x
            ## FIXME: THIS ASSUMES X CONTAINS ALL STATES (don't have age cats otherwise)
            state <- setNames(numeric(length(state_names)), state_names)
        }
    }
    if (vaxify) {
        if (!is.null(params)) {
            state <- expand_state_vax(
                state,
                get_vax_model_type(get_vax(params))
            )
        }
        if (!is.null(x)) {
            state_names <- names(x) ## update state names based on those provided in x
            ## FIXME: THIS ASSUMES X CONTAINS ALL STATES (don't have vax cats otherwise)
            state <- setNames(numeric(length(state_names)), state_names)
        }
    }
    if (testify) state <- expand_stateval_testing(state, method = "untested")


    ## if state vector, x, is not provided, either use given N & E0, or
    ## use eigenvector approach
    if (is.null(x)) {
        if (!use_eigvec) {
            ## just use N & E0

            ## get indexing vector to pull out E states, depending on strata
            E_regex <- case_when(
                vaxify ~ "^E*?_unvax*?", ## works for vaxify, whether or not ageify
                testify ~ "^E_u", ## first testify compartment only
                TRUE ~ "^E" ## works for base, ageify
            )

            ## make boolean subsetting vector based on regex over state names
            E_index <- grepl(E_regex, names(state))

            ## get values to assign to E states, depending on strata
            if (ageify) {
                ## one E value per age group, weighted by age-specific pop-size
                E_values <- smart_round(E0 * params[["Ndist"]])
            } else {
                ## otherwise, just one value, as given
                E_values <- E0
            }

            ## assign E values to E indices in state vector
            state[E_index] <- E_values
        } else {
            ## distribute 'E0' value based on dominant eigenvector
            ## here E0 is effectively "number NOT susceptible"
            ## (repair_names_age does nothing to a non-ageified vector)
            if (vaxify) {
                evec <- get_evec(params,
                    testify = testify,
                    # FIXME: expose this with default TRUE
                    do_hazard = getOption("MP_vax_make_state_with_hazard")
                )
            } else {
                evec <- get_evec(params, testify = testify)
            }
            E_values <- repair_names_age(
                smart_round(evec * E0)
            )
            if (any(is.na(E_values))) {
                state[] <- NA
                return(state)
            }
            if (all(E_values == 0)) {
                if (testify) stop("this case isn't handled for testify")
                E_values[["E"]] <- 1
                warning("initial values too small for rounding")
            }
            state[names(E_values)] <- unname(E_values)
        }
        ## update S classes
        ## make sure to conserve N by subtracting starting number infected

        ## get indexing vector to pull out S states, depending on strata
        S_regex <- case_when(
            vaxify ~ "^S.*_unvax.*", ## works for vaxify, whether or not ageify
            testify ~ "^S_[un]",
            TRUE ~ "^S" ## works for base, ageify
        )

        ## make boolean subsetting vector based on regex over state names
        S_index <- grepl(S_regex, names(state))

        ## get S values, depending on case
        if (testify) {
            if (ageify || vaxify) stop("ageify and/or vaxify don't yet work with testify")
            ## FIXME for testify:  (1) make sure get_evec() actually returns appropriate ratios for S
            ##  class; (2) distribute (N-istart) across the S classes, then round
            ## if A = testing rate and B = test-return rate then
            ##  du/dt = -A*u + B*n  [where u is untested, n is negative-waiting]
            ##        = -A*u + B*(1-u)  [assuming we're working with proportions]
            ##     -> -u*(A+B) +B =0 -> u = B/(A+B)
            ## FIXME: get_evec() should work for S!
            ufrac <- with(as.list(params), omega / ((testing_intensity * W_asymp) + omega))
            S_values <- state_round((N - sum(E_values)) * c(ufrac, 1 - ufrac))
        } else if (ageify || vaxify) {
            ## aggregate over states
            istart <- condense_state(E_values)
            if (vaxify) {
                ## aggregate over vax strata
                if (ageify) {
                    ## by age
                    istart <- condense_vax(istart)
                } else {
                    istart <- rowSums(istart)
                }
            }
            ## just get bare values
            istart <- unname(unlist(istart))
            if (length(params[["N"]]) != length(istart)) stop("length of population size vector in params list (params[['N']]) must either be 1 or number of age categories")
            S_values <- params[["N"]] - istart
        } else {
            S_values <- N - sum(E_values)
        }

        state[S_index] <- S_values
    } else {
        if (length(names(x)) == 0) {
            stop("provided state vector must be named")
        }
        if (length(extra_names <- setdiff(names(x), state_names)) > 0) {
            warning("extra state variables (ignored): ", paste(extra_names, collapse = ", "))
        }
        state[names(x)] <- x
    }

    if (normalize) state <- state / sum(state)

    untestify_state <- state ## FIXME: what is this for??
    attr(state, "epi_cat") <- epi_cat
    class(state) <- "state_pansim"

    ## Give a warning if not all state variables are capital letters
    if (!all(first_letter_cap(names(state)))) {
        warning(
            "Not all state variables are capital letters; ",
            "this can result in failure to correctly check ",
            "if a state variable is negative."
        )
    }

    return(state)
}

##' gradient function for ODE runs
##' @param t time vector
##' @param y state vector
##' @param parms parameter vector
##' @param M rate matrix
gradfun <- function(t, y, parms, M) {
    M <- update_ratemat(M, y, parms)
    foi <- update_foi(y, parms, make_betavec(state = y, parms))
    flows <- col_multiply(M, y) ## faster than sweep(M, y, MARGIN=1, FUN="*")

    g <- colSums(flows) - rowSums(flows)
    return(list(g, foi = foi))
}

##' Run simulation across a range of times
##' @inheritParams do_step
##' @param nt number of steps to take
##' @param ratemat_args additional arguments to pass to \code{make_ratemat}
##' @param step_args additional arguments to pass to \code{do_step}
##' @param M rate matrix
##' @param use_ode solve as differential equation?
##' @param ode_args additional arguments to ode()
##' @importFrom stats rnbinom
##' @examples
##' params <- read_params("ICU1.csv")
##' r1 <- run_sim_range(params)
##' r2 <- run_sim_range(params,use_ode=TRUE)
##' matplot(r1[,"t"],r1[,-1],type="l",lty=1,log="y")
##' matlines(r2[,"t"],r2[,-1],lty=2)
##' @importFrom dplyr left_join
##' @export
run_sim_range <- function(params,
                          state = make_state(params[["N"]], params[["E0"]]),
                          nt = 100,
                          dt = 1,
                          M = NULL,
                          ratemat_args = NULL,
                          step_args = NULL,
                          use_ode = FALSE,
                          ode_args = list()) {
    ## cat("beta0",params[["beta0"]],"\n")
    if (is.null(M)) {
        M <- do.call(make_ratemat, c(list(state = state, params = params), ratemat_args))
    }
    if (use_ode) {
        res <- do.call(
            deSolve::ode,
            c(
                nlist(
                    y = state,
                    times = seq(nt) * dt,
                    func = gradfun,
                    parms = params,
                    M
                ),
                ode_args
            )
        )
        res <- dfs(res)
        if (nrow(res) < nt) {
            time_df <- data.frame(time = 1:nt)
            res <- dplyr::left_join(time_df, res)
        }
        names(res)[1] <- "t" ## ode() uses "time"
    } else {
        ## set up outputs

        ## FOIs
        col_suffix <- if (has_age(params) && has_vax(params)) {
            ## both age and vax
            paste0("_", expand_names(get_age(params), get_vax(params)))
        } else if (has_age(params)) {
            ## just age
            paste0("_", get_age(params))
        } else if (has_vax(params)) {
            ## just vax
            paste0("_", get_vax(params))
        } else {
            ""
        }

        colnames <- paste0("foi", col_suffix)
        ncol <- length(col_suffix)
        foi <- matrix(NA,
            nrow = nt, ncol = ncol,
            dimnames = list(NULL, colnames)
        )

        ## RESULTS ARRAY
        res <- matrix(NA,
            nrow = nt, ncol = length(colnames(M)),
            dimnames = list(
                time = seq(nt),
                state = colnames(M)
            )
        )
        ## initialization
        res[1, names(state)] <- state
        foi[1, ] <- update_foi(state, params, make_beta(state, params))

        ## loop over time steps
        if (nt > 1) {
            for (i in 2:nt) {
                ## update state
                state <- do.call(
                    do_step,
                    c(
                        nlist(state,
                            params,
                            ratemat = M,
                            dt
                        ),
                        step_args
                    )
                )
                ## update foi
                foi[i, ] <- update_foi(state, params, make_beta(state, params))

                if (!identical(colnames(res), names(state))) browser()
                res[i, ] <- state
            }
        }
        res <- dfs(t = seq(nt), res, foi)
    }
    ## need to know true state - for cases with obs error
    attr(res, "state") <- state

    if (any(!is.finite(state))) {
        warning("non-finite values in state, in run_sim_range")
    } else {
        bad_states <- state[state < -sqrt(.Machine$double.eps)]
        if (length(bad_states) > 0) {
            warning(
                "negative state variables (below tolerance) in run_sim_range:",
                paste(sprintf("%s=%1.3g", names(bad_states), bad_states), collapse = "; ")
            )
        } else if (any(state < 0)) {
            warning(
                "End of run_sim_range check: One or more state variables is negative, between -sqrt(.Machine$double.eps) and 0 \n Check following message for details \n ",
                paste(utils::capture.output(print(state)), collapse = "\n")
            )
        }
    }

    return(res)
}

##' construct a Gamma-distributed delay kernel
##' @param prop area under the curve (proportion reported)
##' @param delay_mean mean value
##' @param delay_cv coeff of var
##' @param max_len maximum delay
##' @param tail_crit criterion for selecting maximum delay
##' @importFrom stats pgamma qgamma
## mean = a*s
## sd = sqrt(a)*s
## cv = 1/sqrt(a)
## s = mean/cv^2
## a = 1/cv^2
make_delay_kernel <- function(prop, delay_mean, delay_cv, max_len = ceiling(tail_val), tail_crit = 0.95) {
    gamma_shape <- 1 / delay_cv^2
    gamma_scale <- delay_mean / gamma_shape
    tail_val <- qgamma(tail_crit, shape = gamma_shape, scale = gamma_scale)
    if (max_len < tail_val) {
        warning(sprintf(
            "max_len (%d) is less than qgamma(%f, %1.1f, %1.1f)=%1.1f",
            max_len, tail_crit, gamma_shape, gamma_scale, tail_val
        ))
    }
    pp <- diff(pgamma(seq(max_len + 1), shape = gamma_shape, scale = gamma_scale))
    pp <- pp / sum(pp) ## normalize to 1
    v <- prop * pp
    return(v)
}
