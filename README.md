# IntervalLinearAlgebra.jl: Linear algebra done rigorously

## Lightning talk at JuliaCon 2021

These are the slides for the [talk](https://pretalx.com/juliacon2021/talk/WA7BP8/) "IntervalLinearAlgebra.jl: linear algebra done rigorously" at JuliaCon 2021.

<details>
<summary>Abstract</summary>
Solving linear systems is central in most computational domains, from mathematics to engineering applications. This talk will introduce IntervalLinearAlgebra.jl: a package written in Julia to solve linear systems, with interval or real coefficients, rigorously. That is, producing a set guaranteed to contain the true solution of the original problem. This can be applied to solve problems involving uncertainty propagation or perform self-validated computations.
</details>


## Schedule

The talk will be given on **Friday 30th July 2021 at 19:40 UTC** [convert to your timezone](https://arewemeetingyet.com/2021-07-30/19:40)

To access the event, please [register to JuliaCon](https://juliacon.org/2021/tickets/), it's free!

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

- Start the notebook
  ```julia
  julia> using Pluto; Pluto.run(notebook="ILAjuliacon2021.jl")
  ```

## License

The code is released under MIT license and the text under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) license. Copyright Luca Ferranti 2021.
