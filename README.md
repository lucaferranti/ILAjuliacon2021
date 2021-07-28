# IntervalLinearAlgebra.jl: Linear algebra done rigorously

These are the slides for my [talk](https://pretalx.com/juliacon2021/talk/WA7BP8/) at JuliaCon 2021. 

The talk introduces [IntervalLinearAlgebra.jl](https://github.com/lucaferranti/IntervalLinearAlgebra.jl), a package to perform numerical linear algebra rigorously using interval arithmetic.

If you liked the talk, make sure to check the package itself!

## Talk at JuliaCon 2021

The talk will be given on **Friday 30th July 2021 at 19:40 UTC** [convert to your timezone](https://arewemeetingyet.com/2021-07-30/19:40)

To access the event, please [register to JuliaCon](https://juliacon.org/2021/tickets/), it's free!

<details>
<summary>Abstract</summary>
Solving linear systems is central in most computational domains, from mathematics to engineering applications. This talk will introduce IntervalLinearAlgebra.jl: a package written in Julia to solve linear systems, with interval or real coefficients, rigorously. That is, producing a set guaranteed to contain the true solution of the original problem. This can be applied to solve problems involving uncertainty propagation or perform self-validated computations.
</details>


## Setup

- If you haven't already, [install julia](https://julialang.org/downloads/)
- Open the terminal and navigate to the folder where you want to save the slides
- Clone the repository and enter it
  ```
  git clone https://github.com/lucaferranti/ILAjuliacon2021.git
  ```
  ```
  cd ILAjuliacon2021
  ```
- Start a Julia session from the repository and activate the environment
  ```
  julia --project
  ```
- Instantiate the environment, this will download all missing dependencies if necessary.
  ```julia
  using Pkg; Pkg.instantiate()
  ```
- Start the notebook
  ```julia
  using Pluto; Pluto.run(notebook="ILAjuliacon2021.jl")
  ```

## License

The code is released under MIT license and the text under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) license. Copyright Luca Ferranti 2021.
