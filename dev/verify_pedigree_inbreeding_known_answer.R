source("R/pedigree.R")

# Known-answer case: 1,2 unrelated founders. 3,4 = full sibs (both offspring of 1x2).
# 5 = offspring of a full-sib mating (3x4). Classic expected F(5) = 0.25.
ped_df <- data.frame(
  id     = c(1, 2, 3, 4, 5),
  mother = c(0, 0, 1, 1, 3),
  father = c(0, 0, 2, 2, 4)
)

result <- pedigree_inbreeding(ped_df)
print(result)
cat("\nF(5) expected = 0.25, got =", result$f_pedigree[result$id == 5], "\n")
cat("F(1..4) expected = 0 (founders/non-inbred), got =", result$f_pedigree[result$id %in% 1:4], "\n")
