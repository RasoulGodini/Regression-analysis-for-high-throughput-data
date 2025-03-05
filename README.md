The code can help you loop a regression model over all the columns of a dataframe. It can be used for transcriptome, metabolome, lipidome, and proteome data.
This code uses a log10 transformation. The percent change is calculated based on a base 10 power of the estimation value.
Shapiro-Wilk test is used to assess the normality of each column.
VIF (Variance Inflation Factor) is calculated for all the confounders.

Note: Random data has been applied in this example.
Note: Change the confounders according to your project.
Note: Use your specific linear model (lm, glm, etc.).

I hope this helps, and enjoy!
Rasoul Godini
