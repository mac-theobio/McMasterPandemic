## AGE CATEGORY UTILITIES

##' Construct vector of age category labels
##'
##' @param min minimum age
##' @param max maximum age
##' @param da age bin width
##' @family classic_macpan
##' @export
mk_agecats <- function(min = 0, max = 90, da = 10) {
    s1 <- seq(min, max, by = da)
    ## if we have one-year age groups, don't make hyphenated
    ## categories
    if (da == 1) {
        out <- c(
            as.character(s1[-length(s1)]),
            paste0(s1[length(s1)], "+")
        )
    } else {
        out <- c(
            sprintf("%d-%d", s1[-length(s1)], (s1[-1] - 1)),
            paste0(s1[length(s1)], "+")
        )
    }
    return(out)
}

##' Repair ageified names
##'
##' sometimes age categories get changed from 0-10 to 0.10 and 60+ to 60.
##' replace each substituted period with the correct character
##'
##' @param x a named list or data.frame
##'
##' @return a named list or data.frame
##' @family classic_macpan
##' @export
##'
##' @examples
##' state <- c(S_0.10 = 100, E_60. = 1)
##' repair_names_age(state)
repair_names_age <- function(x) {
    the_names <- names(x)
    ## remove "date" from the names since this isn't a state variable
    the_names <- the_names[!grepl("date", the_names)]
    ## if there aren't any x names that have an age groups, just return x
    if (!any(grepl("_\\d+", the_names))) {
        return(x)
    }

    ## first replace all periods with a hyphen
    names(x) <- sub("\\.", "-", names(x))
    ## then replace hyphens at the end of the name or before an underscore with a +
    names(x) <- sub("-$", "\\+", names(x))
    names(x) <- sub("-_", "\\+_", names(x))

    return(x)
}

##' Aggregate age categories into larger age categories
##'
##' A function that takes a vector of ages or age groups (either numeric or
##' factor) and returns assigned age categories based on desired lower bin edges
##' provided by the user.
##'
##' @param age Numeric or factor vector of ages
##' @param bin_edges_lower list of numeric values for the lower edges of desired
##'   age group bins in panel 2. If ages are of type factor, all user-specified
##'   lower bin edges must match lower bin edges in the age data.
##'
##' @family classic_macpan
##' @export
##'
##' @examples aggregate_agecats(seq(0, 90, by = 1), c(25, 45))
aggregate_agecats <- function(age,
                              bin_edges_lower) {

    ## 1. Set up age df and binning column
    ##

    ## Set up a df so we can use dplyr tools to add the age groups
    df <- data.frame(age = age)

    ## By default, use the age column as the column to bin age
    ## groups over
    binning_colname <- "age"

    ## If ages are factors, create new (numeric) column with which
    ## we will do age group binning
    if (is.factor(age) || is.character(age)) {
        df <- (df %>% mutate(lower_bin_edge = as.numeric(
            stringr::str_extract(age, "^[0-9]*")
        )))

        ## Check that lower bin edges, as given by user, can be found in
        ## the data as lower bin edges
        edge_check <- (bin_edges_lower %in% df$lower_bin_edge)
        if (!all(edge_check)) {
            missing_bins <- bin_edges_lower[!edge_check]
            err_msg <- paste0(
                "The following lower bin edges were not found in the original data: ",
                paste(as.character(missing_bins),
                    collapse = ", "
                ), ". Values provided as lower_bin_edges must all be lower bin edges in the data."
            )
            stop(err_msg)
        }

        ## Changing binning column to lower_bin_edge for
        ## factor data, since it is numeric
        binning_colname <- "lower_bin_edge"
    }

    ## 2. Format bin edges
    ##

    ## Sort given bin edges
    bin_edges_lower <- c(sort(bin_edges_lower))

    ## If zero wasn't given as a lower bin edge, add it
    if (bin_edges_lower[1] != 0) {
        bin_edges_lower <- c(0, bin_edges_lower)
    }

    ## If Inf wasn't given as as upper limit,
    ## Add Inf to bin edges to properly bin data
    bin_edges <- c(bin_edges_lower, Inf)

    ## 3. Make age group labels
    ##

    ## Make human-readible age group labels from lower bin edges
    labels <- c()
    for (i in 1:length(bin_edges_lower)) {
        ## for all but the last age group, make labels in format of
        ## bottom edge age - (next bottom edge age -1)
        if (i < length(bin_edges_lower)) {
            new_label <- paste0(
                bin_edges_lower[i], "-",
                bin_edges_lower[i + 1] - 1
            )
            labels <- c(labels, new_label)
        } else {
            ## for last label, format is "age+"
            new_label <- paste0(bin_edges_lower[i], "+")
            labels <- c(labels, new_label)
        }
    }

    ## Make lookup table for age group labels
    labels_lookup <- data.frame(
        id = 1:length(labels),
        age_group = forcats::as_factor(labels)
    )

    ## 4. Generate vector of age groups
    ##

    df %>%
        ## identify each age (or lower bin edge in the case of factor
        ## ages) with the age group bin it belongs in
        mutate(id = .bincode(get(binning_colname), bin_edges,
            right = FALSE,
            include.lowest = TRUE
        )) %>%
        ## left join lookup table to add our age group labels
        left_join(labels_lookup, by = "id") %>%
        ## pull just the age group col
        pull(age_group) -> age_group

    return(age_group)
}

## STATE UTILITIES

##' expand state vector by age classes and population distribution
##'
##' epidemiological state varies fast, age category varies slowly
##' @param x state vector
##' @param age_cat vector of age categories
##' @param Nvec population distribution (as counts)
##' @examples
##' params <- read_params("PHAC_testify.csv")
##' ss <- make_state(params=params)
##' # ss2 <- expand_state_age(ss)
##' @importFrom purrr map_dfr
##' @importFrom dplyr across
##' @importFrom dplyr relocate
##' @importFrom dplyr everything
##' @family classic_macpan
##' @export
expand_state_age <- function(x, age_cat = mk_agecats(),
                             Nvec = NULL) {
    ## workaround for rcmdcheck
    total <- S <- NULL
    ## save original attributes
    original_attributes <- attributes(x)

    ## if no population is provided, assume a uniform distribution
    if (is.null(Nvec)) {
        Nvec <- mk_Nvec(age_cat, Ntot = sum(x))
    }

    ## expand state labels with ages
    new_names <- expand_names(names(x), age_cat)

    if (sum(Nvec) > 1) {
        ## if we have population vector with counts, go through the
        ## process of distributing counts carefully across subcategories

        ## check that number of age groups matches
        if (length(age_cat) != length(Nvec)) {
            stop("population distribution must have same length as the number of age categories")
        }

        ## check that state vector and population distribution
        ## sum to the same value
        if (sum(x) != sum(Nvec)) {
            stop("state and population distributions must have same sum")
        }

        ## make population distribution
        Ndist <- Nvec / sum(Nvec)

        ## distribute total state counts over subcategories in a way that respects
        ## the overall population distribution
        x <- (
            ## get all states but S
            x[!grepl("^S", names(x))]
            ## distribute counts across all states but S
            ## guided by population distribution
                %>%
                map_dfr(distribute_counts, dist = Ndist)
                ## get total counts by age over all states but S
                %>%
                mutate(total = rowSums(across(everything())))
                ## calculate S by age as the age-specific population - count over all
                ## other states
                %>%
                mutate(S = Nvec - total)
                ## drop total col and move S to far left of the table (so state columns
                ## are ordered as we usually order them)
                %>%
                select(-total)
                %>%
                relocate(S)
        )

        ## flatten to vector, AND ensure counts
        ## are first sorted by age category, then by state category
        x <- as.vector(t(x))

        ## add expanded names
        names(x) <- new_names
    } else {
        x_base <- x ## save old (unexpanded) states before converting x to output
        ## initialize zero state vec (and return this unless we're dealing with a
        ## normalized pop where the total population is 1)
        x <- rep(0, length(new_names))
        ## add expanded names
        names(x) <- new_names
        if (sum(Nvec) == 1) {
            ## if we have a population vector normalized to 1, distribute state
            ## among subcategories according to population distribution
            for (state in names(x_base)) {
                vals <- x_base[[state]] * Nvec
                x[grepl(paste0("^", state), names(x))] <- vals
            }
        }
    }

    ## update names in original attributes to new names
    original_attributes$names <- names(x)

    ## restore original attributes
    attributes(x) <- original_attributes

    ## add sensible names and age cat attribute to output
    attr(x, "age_cat") <- age_cat

    return(x)
}

##' collapse age/testing/etc. categories (only; don't do other condensation)
##' @param x a state vector or simulation result df (an object of class `state_pansim` or `pansim`)
##' @importFrom tidyr separate
##' @family classic_macpan
##' @export
## wrap into condense() ?
## TODO: rewrite this as a generic funtion with custom methods for state_pansim
## and pansim objects
condense_age <- function(x) {
    ## R CMD CHECK doesn't understand dplyr, this is a workaround
    state <- value <- NULL
    ## check if there are even age categories in the object, and if not, return x
    if (!has_age(x)) {
        return(x)
    }

    ## save original attributes
    original_attributes <- attributes(x)

    ## get input type
    input_class <- class(x)

    ## if x is a vector (e.g. a state vec), convert it to a (wide) df
    if (is.null(nrow(x))) {
        x <- as.data.frame(t(unclass(x)))
    }

    ## pivot longer to make state aggregation easier
    ## take care to preserve date column, if it exists
    if ("pansim" %in% input_class) {
        (x
            %>% pivot_longer(-date,
                names_to = "var"
            )
        ) -> x
    } else {
        (x
            %>% pivot_longer(everything(),
                names_to = "var"
            )
        ) -> x
    }

    (x
    ## separate state from subcategory
        %>% separate(var,
            into = c("state", "subcat"),
            sep = "_",
            extra = "merge"
        )
        ## turn state column into factor to keep original ordering
        %>% mutate(state = factor(state, levels = unique(state)))
    ) -> x

    ## aggregate value by state (and timestep, if it exists)
    if ("pansim" %in% input_class) {
        (x
            %>% group_by(date, state)
        ) -> x
    } else {
        (x
            %>% group_by(state)
        ) -> x
    }

    (x
        %>% summarise(value = sum(value), .groups = "drop")
        ## put data back into wide format
        %>% pivot_wider(
            names_from = "state",
            values_from = "value"
        )
    ) -> x

    ## finalize output type
    if ("state_pansim" %in% input_class) {
        x <- unlist(x, use.names = TRUE)
    }

    if ("pansim" %in% input_class) {
        x <- as.data.frame(x)
    }

    ## add back original attributes
    new_names <- names(x)
    original_attributes$names <- new_names
    attributes(x) <- original_attributes

    ## remove age attribute
    attr(x, "age_cat") <- NULL

    return(x)
}

## POPULATION DISTRIBUTION

##' generate a population distribution (as a vector of counts)
##' @param age_cat age categories (made with `mk_agecats()`)
##' @param Ntot total population size
##' @param dist (character) either "unif" for uniform or "rand" for a random population distribution
##' @param names (logical) return a named vector? (names are given age categories)
##' @importFrom stats runif
##' @family classic_macpan
##' @export
mk_Nvec <- function(age_cat = mk_agecats(),
                    Ntot = 1e6,
                    dist = c("unif", "rand"),
                    names = FALSE) {
    dist <- match.arg(dist)
    if (is.null(dist)) dist <- "unif"

    n <- length(age_cat)

    ## uniform distribution
    if (dist == "unif") {
        Ndist <- rep(1 / n, n)
    }

    ## random distribution
    if (dist == "rand") {
        Ndraw <- runif(n)
        Ndist <- Ndraw / sum(Ndraw) ## normalize
    }

    ## if total population size is 1 (normalized pop)
    if (Ntot == 1) {
        ## just return the distribution
        Nvec <- Ndist
    } else {
        ## otherwise, distribute counts carefully, avoiding changes in pop size due
        ## to rounding
        Nvec <- distribute_counts(total = Ntot, dist = Ndist)
    }

    ## add names?
    if (names) {
        names(Nvec) <- age_cat
    }

    return(Nvec)
}

## CONTACT MATRIX

##' generate a contact matrix (where each row is a probability distribution)
##' @param age_cat age categories (made with `mk_agecats()`)
##' @param dist (character) either "unif" for uniform, "rand" for a random distributions by age, or "diag" for strictly within-age mixing (pmat = the identity matrix)
##' @param Nvec population distribution vector (needed to generate a "balanced" random contact probability matrix if population distribution isn't perfectly uniform)
##' @family classic_macpan
##' @export
mk_pmat <- function(age_cat = mk_agecats(),
                    dist = "unif",
                    Nvec = NULL) {
    if (!(dist %in% c("unif", "diag"))) stop("dist must be either 'unif' or 'diag'")
    ## get nrow and ncol for pmat
    n <- length(age_cat)

    ## generate contact abundance by group (and set up correct dimnames)
    ## this needs to be a symmetric matrix
    if (dist == "unif") {
        pmat <- matrix(1, nrow = n, ncol = n, dimnames = list(age_cat, age_cat))
    }

    if (dist == "diag") {
        pmat <- diag(1, nrow = n, ncol = n, names = TRUE)
        dimnames(pmat) <- list(age_cat, age_cat)
    }

    ## row normalize to finish
    pmat <- pmat / rowSums(pmat)

    return(pmat)
}

##' Edit parameter description attribute to include age-specific definitions
##' (helper function for `expand_params_age()`)
##'
##' @param params_desc parameter descriptions, as initialized in `read_params()`
##'
##' @family classic_macpan
##' @export
##'
##' @examples
##' params <- read_params("PHAC_testify.csv")
##' params_desc <- attr(params, "description")
##' expand_params_desc_age(params_desc)
expand_params_desc_age <- function(params_desc) {
    params_desc[["beta0"]] <- paste(params_desc[["beta0"]], "per age category (if single value, this same value is taken for all age groups)")
    params_desc[["N"]] <- paste(params_desc[["N"]], "per age category")
    params_desc[["pmat"]] <- "Contact probability matrix: the entry in row i and column j represents the probability that an individual of age i contacts an individual of age j (but not necessarily the other way around---only matrix rows must sum to 1)"

    return(params_desc)
}

##' Expand parameter list to include age structure
##'
##' @param params parameter list (e.g. read in with `read_params()`)
##' @param age_cat vector of age categories
##' @param beta0 vector of age-specific beta0 values (if NULL, assume same beta0 for all age groups as already provided in params)
##' @param transmissibility proportion of contacts between S & I that lead to transmission
##' @param contact_rate_age average overall contact rate by age
##' @param Nvec population distribution (as counts); default is uniform
##' @param pmat contact matrix; default is uniform
##' @param balance_warning should a warning about the balance of contacts be provided?
##' @examples
##' params <- read_params("PHAC_testify.csv")
##' params_age <- expand_params_age(params)
##' @family classic_macpan
##' @export
expand_params_age <- function(params,
                              age_cat = mk_agecats(),
                              beta0 = NULL,
                              transmissibility = NULL,
                              contact_rate_age = NULL,
                              Nvec = NULL,
                              pmat = NULL,
                              balance_warning = TRUE) {
    ## check inputs
    ##

    ## if all three of beta0, transmissibility, and contact_rate_age are provided
    if (all(!is.null(beta0), !is.null(transmissibility), !is.null(contact_rate_age))) {
        ## check lengths
        if (length(beta0) != length(contact_rate_age)) stop("beta0 and contact_age_rate must have same length")
        ## check that these params are consistently defined
        if (!isTRUE(all.equal(beta0, transmissibility * contact_rate_age))) stop("beta0, transmissibility, and contact_rate_age do not satisfy beta0 = transmissibility*contact_rate_age as provided")
    } else if (xor(is.null(transmissibility), is.null(contact_rate_age))) {
        ## if only one of transmissibility or contact_rate_age is provided
        stop("either both transmissibility and contact_rate_age must be provided, or neither")
    }

    if (!is.null(beta0)) {
        if (!(length(beta0) %in% c(1, length(age_cat)))) stop("beta0 must be either a scalar or a vector length age_cat")
    }

    if (!is.null(transmissibility)) {
        if (!(length(transmissibility) %in% c(1, length(age_cat)))) stop("transmissibility must be either a scalar or a vector length age_cat")
    }

    if (!is.null(contact_rate_age)) {
        if (length(contact_rate_age) != length(age_cat)) stop("contact_rate_age must be a vector length age_cat")
    }

    ## prep inputs
    ##

    ## grab existing parameter descriptions
    params_desc <- attr(params, "description")

    ## convert to list
    params <- as.list(params)

    ## perform updates
    ##

    if (is.null(Nvec)) {
        Nvec <- mk_Nvec(age_cat = age_cat, Ntot = params[["N"]])
    }

    ## set defaults for pmat and Nvec
    if (is.null(pmat)) {
        pmat <- mk_pmat(age_cat = age_cat, Nvec = Nvec)
    }

    ## update pop (counts and distribution)
    params[["N"]] <- Nvec
    params[["Ndist"]] <- Nvec / sum(Nvec)
    ## update pmat
    params <- c(params, list(pmat = pmat))

    ## calculate transmission-related params, depending on what is provided

    ## if both transmissibility and contact_rate_age are provided, update beta0
    if (all(!is.null(transmissibility), !is.null(contact_rate_age))) {
        ## (this is fine to do even if beta0 is provided separately thanks to the
        ## initial check that, when all three parameters are provided, their
        ## definition is consistent)
        names(contact_rate_age) <- age_cat
        beta0 <- transmissibility * contact_rate_age
    }

    ## prep outputs
    ##

    ## attach transmission-based params (may be NULL)
    if (!is.null(beta0)) {
        params[["beta0"]] <- beta0
    }
    params[["transmissibility"]] <- transmissibility
    params[["contact_rate_age"]] <- contact_rate_age

    ## check outputs
    ##

    ## check for balance in contacts (assuming constant transmissibility across
    ## age groups; otherwise, would need to pre-multiply by vector c_i,
    ## age-specific contact rates, instead of beta0_i's)
    if (balance_warning) {
        if (!isSymmetric((params[["beta0"]] * params[["N"]]) * params[["pmat"]])) warning("implied total contact rate between age groups is not balanced for this choice of beta0s, population distribution, and contact matrix (use at your own risk!)")
    }

    ## add attributes
    attr(params, "description") <- expand_params_desc_age(params_desc)
    attr(params, "age_cat") <- age_cat

    ## add class
    class(params) <- "params_pansim"

    return(params)
}

## UTILITIES TO USE MISTRY CONTACT MATRICES

##' Check compatibility of user-specified age categories with age categories in
##' Mistry data
##'
##' @param age_cat age categories (e.g. generated by `mk_agecats()`)
##' @family classic_macpan
##' @export
##'
##' @examples check_age_cat_compatibility(mk_agecats(min = 0, max = 80, da = 2))
check_age_cat_compatibility <- function(age_cat) {
    ## process age cat string into integers
    ages <- as.integer(sub("\\+", "", unlist(strsplit(age_cat, "-")), "+", ""))
    if (min(ages) != 0) stop("minimum age must be 0")
    if (max(ages) > 84) stop("maximum age can't exceed 84")
}

## dput(unique(substr(list.files(system.file("params", "mistry-cmats", package="McMasterPandemic")), 1, 2)))
provterr_abbrevs <- c(
    "AB", "BC", "CA", "MB", "NB", "NL", "NS", "NT", "NU", "ON",
    "PE", "QC", "SK", "YT"
)

##' Make population distribution using Mistry et al. data
##'
##' @inheritParams expand_params_mistry
##'
##' @return a vector of population counts (one per age category)
##' @importFrom readr read_csv cols col_double
##' @family classic_macpan
##' @export
##'
##' @examples mk_mistry_Nvec()
mk_mistry_Nvec <- function(province = "ON",
                           age_cat = NULL) {
    if (nchar(province) > 2) {
        stop("please use two-letter abbreviations for provinces and territories")
    }
    province <- toupper(province)
    if (!province %in% provterr_abbrevs) {
        stop("unknown province territory ", sQuote(province))
    }

    ## workaround for rcmdcheck
    age <- age_group <- pop <- NULL
    ## flag whether any aggregation needs to be done across ages
    aggregate <- !is.null(age_cat)

    ## if age categories are user-specified, check that they specify an
    ## aggregation compatible with Mistry's default age-categories; otherwise, use
    ## the default age cats
    if (aggregate) {
        check_age_cat_compatibility(age_cat)
    } else {
        age_cat <- mk_agecats(min = 0, max = 84, da = 1)
    }

    ## set up filename prefix/suffix for each
    ## setting-specific matrix
    if (is.null(province)) {
        filename_prefix <- paste0("CA_")
    } else {
        filename_prefix <- paste0(
            province,
            "_"
        )
    }

    pop_dist <- (read_csv(system.file("params", "mistry-cmats",
        paste0(
            filename_prefix,
            "age_distribution_85.csv"
        ),
        package = "McMasterPandemic"
    ),
    col_names = c("age", "pop"),
    col_types = cols(
        .default = col_double()
    )
    ))

    if (aggregate) {
        bin_edges_lower <- as.numeric(sub("-\\d+|\\+", "", age_cat))
        pop_dist <- (pop_dist
            %>% mutate(
                age_group
                = aggregate_agecats(age, bin_edges_lower)
            )
            %>% group_by(age_group)
            %>% summarise(pop = sum(pop))
        )
    }

    return(pop_dist$pop)
}

##' Expand parameter list to include age structure using Mistry et al. data
##'
##' @param params baseline parameters loaded from one of the parameter files (??)
##' @param transmissibility probability of transmission upon contact with an infected (beta0 = transmissibility * contact_rate)
##' @param province province for which to construct the matrix (if NULL, make a Canada-wide contact matrix)
##' @param contact_rate_setting named list containing setting-specific contact rates (in units of average contacts in the given setting per individual of age i with individuals of age j per day)
##' @param age_cat (optional) list of age groups to aggregate ages in; use `mk_agecats()` to generate (must start with 0 and end with 84); default is single ages starting with 0 and up to 83, then a single 84+ category
##'
##' @return an object of class `params_pansim`
##' @importFrom dplyr summarise group_by pull
##' @importFrom tidyr pivot_wider pivot_longer
##' @family classic_macpan
##' @export
##'
##' @examples
##' params <- read_params("PHAC_testify.csv")
##' age_cat <- mk_agecats(min = 0, max = 80, da = 10)
##' expand_params_mistry(params = params, province = "ON", age_cat = age_cat)
expand_params_mistry <- function(params,
                                 transmissibility = 1,
                                 province = "ON",
                                 contact_rate_setting = mk_contact_rate_setting(),
                                 age_cat = NULL) {
    ## workaround for rcmdcheck
    row_age_group <- col_age_group <- value <- NULL
    ## CHECK AND PREP INPUTS

    ## flag whether any aggregation needs to be done across ages
    aggregate <- !is.null(age_cat)

    ## if age categories are user-specified, check that they specify an
    ## aggregation compatible with Mistry's default age-categories; otherwise, use
    ## the default age cats
    if (aggregate) {
        check_age_cat_compatibility(age_cat)
    } else {
        age_cat <- mk_agecats(min = 0, max = 84, da = 1)
    }
    ## get number of age categories
    n <- length(age_cat)

    ## check setting-specific contact rate list
    check_contact_rate_setting(contact_rate_setting)

    ## BUILD POPULATION DISTRIBUTION(S)
    Nvec <- mk_mistry_Nvec(province = province, age_cat = age_cat)

    ## if aggregating, set up lower bin edges of age groups,
    ## load population by age, and generate new aggregated population counts
    if (aggregate) {
        ## make Nvec
        Nvec_orig <- mk_mistry_Nvec(
            province = province,
            age_cat = mk_agecats(min = 0, max = 84, da = 1)
        )

        ## get lower age bins for aggregation
        bin_edges_lower <- as.numeric(sub("-\\d+|\\+", "", age_cat))
        Nvec_agg <- Nvec
    }

    ## BUILD CONTACT PROBABILITY MATRIX

    ## combine setting-specific frequency matrices through a
    ## linear combination with the specified contact_rate_setting (in units
    ## of avg number of setting-specific contacts per
    ## individual of age i per unit time) to generate an
    ## overall contact matrix (in units of avg number of
    ## contacts per individual of age i per day)

    ## set up filename prefix/suffix for each
    ## setting-specific matrix
    if (is.null(province)) {
        filename_prefix <- paste0("CA_")
    } else {
        filename_prefix <- paste0(
            province,
            "_"
        )
    }
    filename_suffix <- "_setting_85.csv"

    settings <- c("household", "school", "work", "community")

    ## preallocate memory for the frequency matrix list and the contact
    ## probability matrix
    fmats <- sapply(settings, function(x) NULL)
    pmat <- matrix(rep(0, n * n), nrow = n)
    rownames(pmat) <- age_cat
    colnames(pmat) <- age_cat

    for (set in settings) {
        filename <- system.file("params", "mistry-cmats",
            paste0(
                filename_prefix,
                "F_",
                set,
                filename_suffix
            ),
            package = "McMasterPandemic"
        )

        ## load setting-specific matrix
        set_mat <- readr::read_csv(filename,
            col_names = FALSE,
            col_types = cols(
                .default = col_double()
            )
        )

        ## aggregation step
        ## aggregate matrix entries into age groups (sum blocks)---use a pivot trick?
        ## divide rows by new N_i (aggregated)
        if (aggregate) {

            ## multiply rows by original N_i
            set_mat <- as_tibble(as.matrix(set_mat) * Nvec_orig)

            ## aggregate
            colnames(set_mat) <- mk_agecats(0, 84, 1)

            set_mat <- (set_mat
                %>% mutate(row_age_group = aggregate_agecats(
                    mk_agecats(0, 84, 1),
                    bin_edges_lower
                ))
                %>% pivot_longer(-row_age_group,
                    names_to = "col_age_group"
                )
                %>% mutate(col_age_group = aggregate_agecats(
                    col_age_group,
                    bin_edges_lower
                ))
                %>% group_by(row_age_group, col_age_group)
                %>% summarise(value = sum(value), .groups = "drop")
                %>% pivot_wider(names_from = col_age_group)
                %>% select(-row_age_group)
                %>% as.matrix()
            )

            ## row normalize by aggregated population
            set_mat <- set_mat / Nvec_agg
        }

        ## convert to matrix if need be
        if (!aggregate) set_mat <- as.matrix(set_mat)

        ## save setting frequency matrices in case we want to fiddle with contact_rate_setting
        ## later
        rownames(set_mat) <- age_cat
        colnames(set_mat) <- age_cat
        fmats[[set]] <- set_mat

        ## update overall contact matrix by adding a weighted
        ## version of the current setting-specific frequency
        ## matrix
        pmat <- pmat + contact_rate_setting[[set]] * set_mat
    }

    ## row-normalize overall contact matrix to make it a probability matrix
    ## capture rowsums contact matrix (avg contact rate by age)
    contact_rate_age <- rowSums(pmat)
    pmat <- pmat / contact_rate_age

    ## PREPARE OUTPUT

    ## attach standard components to params list
    params <- expand_params_age(
        params,
        age_cat = age_cat,
        pmat = pmat,
        Nvec = Nvec,
        transmissibility = transmissibility,
        contact_rate_age = contact_rate_age
    )
    ## save parameter descriptions to append after adding to the params list
    params_desc <- attr(params, "description")

    ## attach age- and Mistry-specific components to params list
    params <- c(
        params,
        list(
            mistry_contact_rate_setting = contact_rate_setting,
            mistry_fmats = fmats
        )
    )

    ## attach attributes
    attr(params, "description") <- c(params_desc,
        mistry_contact_rate_setting = "Average number of contacts per setting across all age groups (assumed); calculated using Mistry et al. 2021 contact matrices",
        mitsry_fmats = "Setting-specific contact frequency matrices, where row i gives the contact frequency per susceptible of age group i (assumed); from Mistry et al. 2021"
    )
    attr(params, "age_cat") <- age_cat
    class(params) <- "params_pansim"

    return(params)
}

##' Helper function to generate a vector of setting-specific contact rates
##'
##' @param values_to_update named list of values to update in default list of setting-specific contact rates
##'
##' @return named list of setting-specific contact rates
##' @family classic_macpan
##' @export
##'
##' @examples
##' mk_contact_rate_setting()
##' mk_contact_rate_setting(list(school = 0))
mk_contact_rate_setting <- function(values_to_update = NULL) {
    ## initialize with default setting-specific weights
    contact_rate_setting <- list(
        household = 4.11, school = 11.41,
        work = 8.07, community = 2.79
    )

    if (!is.null(values_to_update)) {
        contact_rate_setting <- update_contact_rate_setting(
            contact_rate_setting,
            values_to_update
        )
    }

    return(contact_rate_setting)
}

##' Helper function to update setting-specific contact rate list
##'
##' @inheritParams expand_params_mistry
##' @inheritParams mk_contact_rate_setting
##'
##' @return named list of setting-specific contact rates
##' @family classic_macpan
##' @export
##'
##' @examples
##' update_contact_rate_setting(mk_contact_rate_setting(), list(school = 0))
update_contact_rate_setting <- function(contact_rate_setting,
                                        values_to_update) {
    ## check initialization
    check_contact_rate_setting(contact_rate_setting)

    ## check update values are properly initialized
    if (!is.list(values_to_update)) stop("update values must be provided as a (named) list")
    if (!all(names(values_to_update)
    %in% names(contact_rate_setting))) {
        stop("names in values to update list must be 'household', 'school', 'work', or 'community'")
    }

    ## perform updates
    for (setting in names(values_to_update)) {
        contact_rate_setting[[setting]] <- values_to_update[[setting]]
    }

    return(contact_rate_setting)
}

##' Helper function to check setting-specific contact rate list initialization
##'
##' @param contact_rate_setting named list containing setting-specific contact rates (in units of average contacts in the given setting per individual of age i with individuals of age j per day)
##'
##' @family classic_macpan
##' @export
##' @examples
##' check_contact_rate_setting(mk_contact_rate_setting())
check_contact_rate_setting <- function(contact_rate_setting) {
    ## check that it's a list
    if (!is.list(contact_rate_setting)) stop("setting-specific contact rate must be specified as a list")
    ## check that the avg_contact_rate_per_setting list is correctly specified
    if (!isTRUE(all.equal(
        sort(names(contact_rate_setting)),
        c("community", "household", "school", "work")
    ))) {
        stop("setting-specific contact rate list must have names 'household', 'school', 'work', 'community'")
    }
}

##' Update Mistry-based parameters
##'
##' Recalculate beta0 and/or pmat if setting-specific contact rates and/or
##' overall transmissiblity are changed
##'
##' @param params parameter list initialized using `expand_params_mistry()`
##' @inheritParams expand_params_mistry
##'
##' @details
##' For `contact_rate_setting`, a full list (initialized with `mk_contact_rate_setting()`) or a partial list (e.g. `list(community = 0)`) can be provided. If a partial list is provided, default values from Mistry et al. 2021 are assumed for the contact rates that have not been specified (see `mk_contact_rate_setting()` for values).
##'
##' @return an object of class `params_pansim`
##' @family classic_macpan
##' @export
##'
##' @examples
##' params <- read_params("PHAC_testify.csv")
##' age_cat <- mk_agecats(min = 0, max = 80, da = 10)
##' params_mistry <- expand_params_mistry(params = params, province = "ON", age_cat = age_cat)
##' params_mistry <- update_params_mistry(params_mistry,
##'  contact_rate_setting = list(school = 0))
update_params_mistry <- function(params,
                                 contact_rate_setting = NULL,
                                 transmissibility = NULL) {
    ## perform updates based on new contact_rate_setting
    if (!is.null(contact_rate_setting)) {
        ## checks
        if (!("mistry_contact_rate_setting" %in% names(params))) stop("parameters must first be initialized by expand_params_mistry()")

        ## update from full or partial list of setting-specific contact rates
        contact_rate_setting <- update_contact_rate_setting(
            params[["mistry_contact_rate_setting"]],
            contact_rate_setting
        )

        ## get frequency matrices
        fmats <- params[["mistry_fmats"]]

        ## update pmat
        age_cat <- attr(params, "age_cat")
        n <- length(age_cat)
        pmat <- matrix(rep(0, n * n),
            nrow = n,
            dimnames = list(age_cat, age_cat)
        )

        for (setting in names(contact_rate_setting)) {
            pmat <- pmat + contact_rate_setting[[setting]] * fmats[[setting]]
        }

        ## get implied age-specific contact rates and row-norm pmat
        contact_rate_age <- rowSums(pmat)
        pmat <- pmat / contact_rate_age

        ## update params entries
        params[["mistry_contact_rate_setting"]] <- contact_rate_setting
        params[["contact_rate_age"]] <- contact_rate_age
        params[["pmat"]] <- pmat
    }

    ## update transmissibility
    if (!is.null(transmissibility)) {
        params[["transmissibility"]] <- transmissibility
    }

    ## regenerate beta0
    params[["beta0"]] <- with(params, transmissibility * contact_rate_age)

    return(params)
}

## SIMULATION

##' Ageify a basic simulation
##'
##' Wrapper function for `run_sim()` that incorporates existing ageify tools
##'
##' @param base_params base parameter set without age structure (e.g. generated by `read_params()`)
##' @param base_state base state vector without age structure (e.g. generated by `make_state()`)
##' @param age_cat age categories (e.g. generated by `mk_agecats()`)
##' @param beta0 vector of beta0 values for each age group
##' @param pmat contact matrix (e.g. generated by `mk_pmat()`)
##' @param Nvec population distribution (e.g. generated by `mk_Nvec()`)
##' @param ... additional arguments to `run_sim()`
##'
##' @return simulation results (an object of class `pansim`)
##' @family classic_macpan
##' @export
##'
##' @examples
##' params <- update(read_params("PHAC_testify.csv"), testing_intensity=0)
##' state <- make_state(params=params)
##' run_sim_ageify(base_params = params, base_state = state)
run_sim_ageify <- function(base_params,
                           base_state,
                           age_cat = mk_agecats(),
                           beta0 = NULL,
                           pmat = mk_pmat(),
                           Nvec = mk_Nvec(),
                           ...) {
    params <- expand_params_age(base_params,
        age_cat = age_cat,
        beta0 = beta0,
        pmat = pmat,
        Nvec = Nvec
    )
    state <- expand_state_age(base_state,
        age_cat = age_cat,
        Nvec = Nvec
    )

    res <- run_sim(params, state, ...)

    return(res)
}

## 3
## ben's old code

##' @rdname expand_state_age
##'
##' @param ratemat rate matrix
##' @param params parameter vector
##' @examples
##' params <- read_params("PHAC_testify.csv")
##' state <- make_state(params=params)
##' M <- make_ratemat(state,params, sparse=TRUE)
##' Ma <- ageify(M, params)
##' library(Matrix)
##' Matrix::image(Ma)
##' Mta <- ageify(testify(M,params),params)
##' Matrix::image(Mta)
##' @family classic_macpan
##' @export
ageify <- function(ratemat, params, age_cat = mk_agecats()) {
    m <- Matrix::kronecker(diag(length(age_cat)), ratemat)
    new_names <- unlist(lapply(rownames(m), paste, age_cat, sep = "_"))
    dimnames(m) <- list(new_names, new_names)
    return(m)
}

## FIXME: age-dependent params??

if (FALSE) {
    devtools::load_all()
    params <- read_params("PHAC_testify.csv")
    ss <- make_state(params = params)
    ss2 <- expand_state_age(ss)
    ## hack so we have an infective
    ss2["Im_11-20"] <- 1
    ss2["E_91+"] <- 0
    tot_I <- function(x) sum(x[grep("^I[a-z]", names(x))])
    tot_I(ss)
    sum(ss2)
    tot_I(ss2)
    condense.pansim(data.frame(date = NA, rbind(ss2)), add_reports = FALSE)
    M <- make_ratemat(ss2, params, sparse = TRUE)
    show_ratemat(M)
    aa <- mk_agecats()
    ## compound symmetric example
    pmat <- matrix(0.1, nrow = length(aa), ncol = length(aa), dimnames = list(aa, aa))
    diag(pmat) <- 1
    ifun <- function(M) {
        Matrix::image(Matrix(M), scales = list(
            y = list(at = seq(nrow(M)), labels = rownames(M)),
            x = list(at = seq(ncol(M)), labels = colnames(M), rot = 90)
        ))
    }
    params_age <- c(as.list(params), list(pmat = pmat))
    b1 <- make_beta(ss2, params_age, full = FALSE)
    ifun(b1)
    b2 <- Matrix(make_beta(ss2, params_age, full = TRUE))
    ifun(b2)
    M <- make_ratemat(ss2, params_age, sparse = TRUE)
    show_ratemat(M)
    M %*% ss2
}
