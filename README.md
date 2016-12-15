# grbmet

GBRMET simulates regulated Brownian Motion and generates escape times through a target boundary. The parameter M in the code represents the location of the reflecting boundary. Such a Brownian Motion is called "regulated" or "reflected". When M=inf, the model is the standard Brownian Motion. When M=finite, Browninan motion will be absorbed at -M when it attempts to go below it. 

The escape times follow a power-law tail in the case of free Brownian Motion (M=inf) and a bounded, intermediate, truncated --whatever you call it-- power-law in the case of Reflected Brownian Motion (M=finite). The power law exponent in each case is 1.5 in the limit.