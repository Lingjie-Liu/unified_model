library(tseries)
?acf()

data <- c(1, -1, 0.5, -0.5, 0.7, -0.8)
data <- ts(data)
plot.ts(data)

print(acf(data))
