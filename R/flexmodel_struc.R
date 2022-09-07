if (FALSE) {
xx = epi_mat_list(
  fred = 1,
  tom = matrix(
    const_named_vector(letters[7:12], 0),
    3, 2,
    dimnames = list(letters[1:3], letters[4:5])
  ),
  sam = matrix(
    letters[1:6],
    3, 2,
    dimnames = list(letters[1:3], letters[4:5])
  )
) %>% derive(dave = ~ sam %*% t(sam), jen = ~ sam * tom, carl = ~ fred)
names(xx$tom)
values(xx$tom)
names(xx$sam)
as.character(xx$sam)
}
