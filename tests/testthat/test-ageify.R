library(McMasterPandemic)
## devtools::load_all()
library(testthat)

context("ageify")

pp <- update(read_params("PHAC_testify.csv"), testing_intensity=0)
ss <- make_state(params=pp)
tot_I <- function(x) sum(x[grep("^I[a-z]",names(x))])

## not really proper tests yet: FIXME/clean me up!
test_that("generic age stuff", {
    ss2 <- expand_stateval_age(ss)
    ## hack so we have an infective
    ss2["Im_11-20"] <- 1
    ss2["E_91+"] <- 0
    expect_equal(sum(ss),sum(ss2))
    expect_equal(tot_I(ss2),1)
    condense.pansim(data.frame(date=NA,rbind(ss2)),add_reports=FALSE)
    M <- make_ratemat(ss2, pp, sparse=TRUE)
    show_ratemat(M)
    aa <- mk_agecats()
    ## compound symmetric example
    Cmat <- matrix(0.1, nrow=length(aa), ncol=length(aa), dimnames=list(aa,aa))
    diag(Cmat) <- 1
    ifun <- function(M) {
        Matrix::image(Matrix(M),scales=list(y=list(at=seq(nrow(M)),labels=rownames(M)),
                                            x=list(at=seq(ncol(M)),labels=colnames(M), rot=90)))
    }
    ppa <- c(as.list(pp),list(Cmat=Cmat))
    b1 <- make_betavec(ss2, ppa, full=FALSE)
    expect_equal(dim(b1), c(10,40))
    ifun(b1)
    b2 <- Matrix(make_betavec(ss2, ppa, full=TRUE))
    ifun(b2)
    expect_equal(dim(b2), c(10,130))
    M <- make_ratemat(ss2, ppa, sparse=TRUE)
    expect_equal(dim(M),c(130,130))
    show_ratemat(M)
    M %*% ss2
    rr <- run_sim_range(ppa, ss2, nt=1000)
    ## suppress warnings about values <0 (??)
    suppressWarnings(matplot(rr[,1],rr[,-1],lty=1,type="l",log="y"))
    rr2 <- run_sim(ppa, ss2,end_date="2022-Nov-1",condense=FALSE)
    plot(rr2,log=TRUE)+theme(legend.position="none")
})
