x <- seq(-10, 10, by = .1)

# Choose the mean as 2.5 and standard deviation as 0.5.
y <- dnorm(x, mean = 2.5, sd = 2)

y1 <- scale(y)
y2 <- (y - min(y)) / (max(y) - min(y))

df <- data.frame(x = x, y = y)
df1 <- data.frame(x = x, y= y1)
df2 <- data.frame(x = x, y= y2)

ggplot()+
  geom_line(data = df, aes(x = x, y = y), color = "blue") +
  geom_line(data = df1, aes(x = x, y = y1), color = "red") +
  geom_line(data = df2, aes(x = x, y = y2), color = "black") 
