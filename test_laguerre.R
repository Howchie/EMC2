library(statmod)

# Gauss-Laguerre nodes and weights
glag <- gauss.quad(20, kind="laguerre")
nodes <- glag$nodes
weights <- glag$weights

print(nodes)
print(weights)
