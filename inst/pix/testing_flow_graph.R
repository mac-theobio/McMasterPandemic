library(diagram)
version <- 1
## version 2 (unfinished) implements DJDE picture with L-to-R flow and no epi transitions

if (version == 1) {
    m <- matrix(
        byrow = TRUE, ncol = 2,
        c(
            0.25, 5 / 6, ## - test
            0.5, 2 / 6, ## untested
            0.25, 3 / 6, ## neg sample
            0.75, 3 / 6, ## pos sample
            0.5, 4 / 6, ## pos test
            0.5, 5 / 6
        )
    ) ## + test
    m[, 2] <- 1 - m[, 2] ## flip y axis
} else {
    m <- matrix(
        byrow = TRUE, ncol = 2,
        c(
            0.5, 1 / 6, ## - test
            0.2, 0.5, ## untested
            0.25, 2 / 6, ## neg sample
            0.25, 4 / 6, ## pos sample
            0.75, 4 / 5, ## pos test
            0.5, 5 / 6
        )
    ) ## + test
}

## dummies appear to be necessary to get the right limits??
mx <- cbind(rep(0:1, each = 3), rep(2:4, 2) / 6)
m <- rbind(m, mx)

labs <- c(
    "negative\nreports",
    "untested",
    "negative\nsample\n(waiting)",
    "positive\nsample\n(waiting)",
    "tested\npositive",
    "positive\nreports"
)

if (version == 2) labs[2] <- c("untested or tested negative")
labs <- c(labs, rep("", 6))
##
A <- matrix(0, 12, 12)
A[3, 2] <- 0.5 ## untested -> neg
A[4, 2] <- 0.5 ## untested -> pos
A[2, 3] <- 1 ## neg back to untested
A[5, 4] <- 1 ## pos to tested positive
## A[cbind(2:5,c(7,8,8,9))] <- 0.9
## A[cbind(c(10,11,11,12),2:5)] <- 0.9

box.type <- c("rect", rep("ellipse", 4), "rect", rep("rect", 6))
box.size <- c(0.075, 0.075, 0.1, 0.1, 0.075, 0.075, rep(0, 6))
box.prop <- rep(0.75, 12)
arr.col <- A
storage.mode(arr.col) <- "character"
arr.col[] <- ifelse(A == 0.9, "blue", "black")

if (version == 2) {
    box.prop[2] <- 0.2
    box.size[2] <- 0.3
}

p <- function(add = FALSE) {
    plotmat(A,
        pos = m, name = labs,
        curve = -0.1, ## flip sign to get neg -> untested above untested -> neg
        arr.col = arr.col,
        arr.lwd = 3,
        box.size = box.size,
        box.type = box.type,
        box.prop = box.prop,
        cex.txt = 0, ## no arrow text
        add = add
    )
}
p()

## arrows for epi transitions
if (version == 1) {
    for (i in c(2, 3, 5)) {
        straightarrow(from = c(0, m[i, 2]), to = m[i, ], lcol = "blue", arr.pos = 0.4)
    }
    straightarrow(from = c(0.55, m[4, 2]), to = m[4, ], lcol = "blue", arr.pos = 0.25)
    for (i in c(2, 4, 5)) {
        straightarrow(from = c(m[i, ]), to = c(1, m[i, 2]), lcol = "darkblue", arr.pos = 0.6)
    }
    straightarrow(from = m[3, ], to = c(0.45, m[3, 2]), lcol = "darkblue", arr.pos = 0.75)
}

curvedarrow(c(0.37, 0.62), m[1, ], lty = 2, curve = 0.5)
curvedarrow(c(0.62, 0.38), m[6, ], lty = 2, curve = -0.3)
p(add = TRUE) ## repeat plot to mask arrow stuff
