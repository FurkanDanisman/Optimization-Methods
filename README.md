# Optimization-Methods

For optimization, by using R, we can have various tools: 

- Newton Raphson Method
- Basic Stochastic Optimization
- Monte Carlo Optimization
- Simulated annealing

None of them are flawless, and each has negative and positive aspects. In this study, I will highlight each method's differences, advantages and disadvantages. However, first, let me introduce the equation to be optimized. 

<p align="center">
  $𝜑(𝑥) = sin(𝜋𝑥^3)𝑒^{−𝑥^2}$          ; $𝑥 ∈ [−2, 2]$  
</p>

The reason why this function was chosen is merely due to its interesting plot and challenging optimization
procedure. To illustrate my reasons, one can observe the below plot with its multiple local maximum
points, which makes it hard to distinguish the global maximum for each method. With its complex shape,
it is hoped that each method's weakness will be highlighted more clearly. Also, it resembles a heart signal
for people like me who have no idea what an actual heart signal looks like.

<p align="center">
  <img width="500" height="250" src="https://github.com/FurkanDanisman/Optimization-Methods/raw/main/images/Untitled.png">
</p>



You may find the analysis and visualization [here.](https://furkandanisman.github.io/Optimization-Methods/Code/Optimization-Methods.html) 
