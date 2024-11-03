# functions.R

### Reference datasets

# Data from the titanic disaster
titanic <- rbind(cbind(rep("Female",  89), rep("Deceased",  89)),
                 cbind(rep("Male",   483), rep("Deceased", 483)),
                 cbind(rep("Female", 230), rep("Survived", 230)),
                 cbind(rep("Male",   112), rep("Survived", 112)))
titanic <- titanic[sample(nrow(titanic)), ]  # Shuffle rows
titanic_df <- data.frame(Sex = titanic[, 1], Condition = titanic[, 2])

# Discretization example from tutorial P03_Funciones
p03_disc_values <- c(11.5, 10.2, 1.2, 0.5, 5.3, 20.5, 8.4)
p03_disc_bins <- 4
cat('p03_disc_values: (', p03_disc_values, '), p03_disc_bins:',
    p03_disc_bins,'\n')

# Data known to have a sample variance of 23.5
sample_variance_23p5_df <- data.frame(data = c(3.0, 4.0, 7.0, 12.0, 14.0))
cat('sample_variance_23p5_df: (',
    paste(sample_variance_23p5_df$data, collapse=', '), ')\n')

# Data known to have a population variance of 2.917
pop_variance_2p917_df <- data.frame(die = c(1, 2, 3, 4, 5, 6))
cat('pop_variance_2p917_df: (',
    paste(pop_variance_2p917_df$die, collapse = ', '), ')\n')

# Data known to give an AUC of 1.0 with ROC
auc_1p0_df <- data.frame(
  Score = c(0.1, 0.4, 0.35, 0.8, 0.5),
  Label = c(FALSE, TRUE, FALSE, TRUE, TRUE)
)

# Data known to give an AUC of 0.75 with ROC
auc_0p75_df <- data.frame(
  Score = c(0.66, 0.09, 0.38, 0.27, 0.81, 0.44, 0.81, 0.81, 0.79, 0.43),
  Label = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE)
)

# Data from tutorial P03, exercise 4, known to have entropy 0.971
p03_entropy_0p971 <- c('a', 'a', 'c', 'c', 'c')
cat('p03_entropy_0p971: (',
    paste(p03_entropy_0p971, collapse = ', '), ')\n')

# Dataframe with mixed column type for column-wise testing
std_mixed_df <- data.frame(
  A = c(1L, 1L, 1L, 3L, 3L),
  B = c(TRUE, TRUE, TRUE, FALSE, FALSE),
  C = c(4.5, 5.5, 6.5, 5.5, 6.5),
  D = c('a', 'a', 'a', 'c', 'c'),
  E = c(0.1, 0.4, 0.35, 0.8, 0.5),
  F = c(FALSE, TRUE, FALSE, TRUE, TRUE),
  G = c(4.5, 5.5, 6.5, 6.5, 6.5),
  H = factor(c('a', 'a', 'a', 'c', 'c'))
)

# Data known to give a correlation of -0.6847868
corr_m0p685_df <- data.frame(
  x = c(   8,    3,    5,    7,    1,    2,    6,    7,    4,    9 ),
  y = c( 2.0,  2.0,  1.5,  1.0,  2.5,  3.0,  1.5,  2.0,  2.0,  1.5 )
)
