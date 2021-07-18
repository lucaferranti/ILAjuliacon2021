### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 0bc12fe2-4d08-4840-b651-b58df23f0277
begin
	using Pkg;
	Pkg.activate(".")
	using IntervalLinearAlgebra, LazySets, Plots, Images, StaticArrays, PlutoUI
end

# ╔═╡ f403b050-e2e6-11eb-1c44-911de0d748bd
html"<button onclick='present()'>present</button>"

# ╔═╡ c9e05bd9-1d7b-41c0-99dc-29b1c1e20090
html"""
<style>
.flex-container {
	display : flex;
}

.left {
	float : left;
	
}

.right {
	float : right;
}
</style>
<h1>IntervalLinearAlgebra.jl</h1>
<h3>Linear algebra done rigorously</h3>
<div class="flex-container">
<div class="left">
<br><br>
<p><b>Luca Ferranti</b>, David Sanders, Marcelo Forets</p>
<p><i>Google Summer of Code (GSoC) 2021</i></p>
<p>Presentation for <i>JuliaCon 2021</i></p></div>
<div class="right">
<img width="450" src="https://raw.githubusercontent.com/lucaferranti/IntervalLinearAlgebra.jl/main/docs/src/assets/logo.png"></div></div>"""

# ╔═╡ c3a5fd72-9a9d-45f3-aaf3-bd61637d286c
md"""
## Interval matrices

An interval matrix $\mathbf{A}\in\mathbb{I}\mathbb{R}^{m\times n}$ is defined as

$\mathbf{A} = \{A \in \mathbb{R}^{m\times n} | A_{ij}\in\mathbf{A}_{ij}\quad i=1,\ldots,m\quad j=1\ldots,n\}$

for example
"""

# ╔═╡ f2cff36e-9136-4eb1-acff-c97af31c2450
AA = [1..2 3..4; 5..6 7..8]

# ╔═╡ 286f40a2-840d-431b-a531-fe9e09dc55f1
md"""
##

 $\mathbf{A}$ can also be represented in *midpoint-radius* notation $A_c\pm A_\Delta$, where $A_c$ is the *midpoint matrix* and $A_\Delta$ is the *radius matrix*.
"""

# ╔═╡ 97e498f5-2e0c-4139-9774-86f58d59d511
AA

# ╔═╡ bb878728-d4ca-4c01-a4f7-bf1811a8f5c2
mid.(AA)

# ╔═╡ a89ee6e9-565a-46d5-82db-bafb251a2200
radius.(AA)

# ╔═╡ 4848561e-2336-40fb-a5c6-8705da85c2ee
md"""
##

### Regular matrices
We say that an interval matrix $\mathbf{A}$ is **regular** if all $A\in\mathbf{A}$ are invertible. Otherwise, the interval matrix is **singular**.

- In general, checking for regularity or singularity is computationally expensive.
"""

# ╔═╡ 5b111f71-4d41-43cf-a3cf-e2a6f4c2adba
html"
<table>
  <tr>
    <th>real</th>
    <th>interval</th>
	<th>check property</th>
  </tr>
  <tr>
    <td>invertible</td>
    <td>regular</td>
	<td>coNP-complete</td>
  </tr>
  <tr>
    <td>singular</td>
    <td>singular</td>
	<td>NP-complete</td>
  </tr>
"

# ╔═╡ 395b75e5-e76f-4f3a-b51f-57eaa3de2d4e
md"""
## Interval linear systems

We are now ready to study the square interval linear system $\mathbf{Ax}=\mathbf{b}$ with $\mathbf{A}\in\mathbb{IR}^{n\times n}$ and $\mathbf{b}\in\mathbb{IR}^n$.
"""

# ╔═╡ f39f94f7-59ff-4dd0-93dd-4f9f158c673c
md"""
The solution set $\mathbf{x}$ is defined as

$\mathbf{x}=\{x\in\mathbb{R}^n | Ax=b\text{ for some } A\in\mathbf{A}, b\in\mathbf{b}\}$
"""

# ╔═╡ 3f5a6fed-6175-4459-8be2-9a330cedc938
md"""
- If $\mathbf{A}$ is regular, then the solution set is non-empty and bounded.
"""

# ╔═╡ c223f95a-9885-48ae-b939-4f3fa54a7987
md"""
##
### Finding the solution set

- How do we solve the square interval linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$?
- A naive approach might be to use Monte Carlo, i.e. generate several random real instances from the interval system.
- We will use the following as running example
"""

# ╔═╡ 32f0ff22-d42f-410f-a1a6-b8ea77b696ff
A = [2..4 -2..1;-1..2 2..4]

# ╔═╡ a79e27c9-399c-44e1-ac79-f65480fe5178
b = [-2..2, -2..2]

# ╔═╡ ba3a1ca6-04dc-45db-b686-a036519f2bc9
md"""
##
"""

# ╔═╡ 38883c77-d506-4e81-9d25-0ac78323241a
begin
	Ns = 100_000
	Astatic = @SMatrix [2..4 -2..1;-1..2 2..4]
	bstatic = @SVector [-2..2, -2..2]
	xs = [rand.(A)\rand.(b) for _ in 1:Ns]
	x = [xs[i][1] for i in 1:Ns]
	y = [xs[i][2] for i in 1:Ns]
	histogram2d(x, y, ratio=1)
end

# ╔═╡ 4549ca50-615c-4fba-afde-303ad231f648
md"""
- Pictures obtained by uniformly samplying 100_000 times each interval.
- **Did we cover the whole set?**
"""

# ╔═╡ 0f196850-f4b1-4830-9198-d72a281da05c
md"""
## Solution set characterization

An important theorem, known as **Oettli-Präger theorem**, helps us describing the solution set $\mathbf{x}$.

###### Oettli-Präger theorem
Given the square interval linear system $\mathbf{Ax}=\mathbf{b}$, then

$
x\in\mathbf{x}\Leftrightarrow|A_cx-b_c|\le A_\Delta|x|+b_\Delta$

- The absolute value is taken elementwise.

- We can remove the absolute values by considering one orthant at the time, obtaining $2^n$ systems of linear inequalities.

- Alternatively, the non-linear inequalities can be solved directly using `IntervalConstraintProgramming.jl`

- Both impractical for higher dimensions!
"""



# ╔═╡ d5e3856d-07e7-4e50-9e4d-c0afb7091d8b
md"""
##
"""

# ╔═╡ 95edbbbf-2e35-446d-9eeb-a53f5c6b10ec
polytopes = solve(A, b, LinearOettliPrager())

# ╔═╡ 4880c8a6-de84-4e55-be90-f73256b0888a
begin
	plot(polytopes, ratio=1)
	histogram2d!(x, y)
end

# ╔═╡ 75a66d2a-01b9-4d6a-8fed-8d8159ba6354
md"""
##
### A 3D example
Below the solution of the interval system with

$
\begin{bmatrix}
[4.5, 4.5]&[0, 2]&[0, 2]\\
[0, 2]&[4.5, 4.5]&[0, 2]\\
[0, 2]&[0, 2]& [4.5, 4.5]
\end{bmatrix}\mathbf{x}=\begin{bmatrix}[-1, 1]\\
[-1, 1]\\
[-1, 1]\end{bmatrix}$

"""

# ╔═╡ cd558a80-2825-4953-870c-93551f60c10e
load("./3dexample.png")

# ╔═╡ c5aa282c-c82b-4521-994d-235cb7934c83
md"""
##
### Exact solution, conclusions
* In general, the solution set of an interval linear system is a non-convex polytope.
  - however it is convex in each orthant.
* The computational complexity to find the solution set grows exponentially with the dimension,
* An alternative more feasible approach is to find a *tight enclosure* of the solution.
"""

# ╔═╡ 32481e35-1a98-4b65-b925-264d51d24743
md"""
## Enclosures of interval linear systems

- We say that an interval box $\mathbf{\Sigma}$ is an enclosure of the set $\mathbf{x}$ if $\mathbf{x}\subseteq\mathbf{\Sigma}$

- Ideally, we want the *hull* of the solution $\mathbf{\Sigma}_{H}$, that is the tightest interval box.

- However, finding the exact hull is in general NP-hard
"""

# ╔═╡ 58b8b1b8-93e7-46c0-b0a0-28532b152908
begin
	plot(interval_hull(ConvexHullArray(polytopes)), label="hull", α=0.3)
	plot!(polytopes, ratio=1, α=1)
end

# ╔═╡ faa03639-a431-4cbf-ac03-bbb33fb7b5b8
md"""
## Algorithms to find the enclosure of the system

`InteralLinearAlgebra.jl` has several algorithms to compute an enclosure of $\mathbf{Ax}=\mathbf{b}$ and an user friendly interface to choose what algorithm and precondition mechanisms to use.

#### Implemented algorithms
- Gaussian elimination
- Gauss-Seidel
- Jacobi
- Hansen-Bliek-Rohn
- Krawczyk
"""

# ╔═╡ 5e0cf15b-21ab-45bb-81fa-ee962a46004f
Xge = solve(A, b, GaussianElimination(), NoPrecondition())

# ╔═╡ fc08d410-02fa-480b-a999-841476dfaaf6
md"""
##
"""

# ╔═╡ 57a5fb2a-5cb2-412f-a1a0-5d27e0a5b8f7
begin
	plot(interval_hull(ConvexHullArray(polytopes)), label="hull", α=0.3)
	plot!(IntervalBox(Xge), label="Gaussian elimination", α=0.2, legend=:right)
	plot!(polytopes, ratio=1, α=1)
end

# ╔═╡ a70a6137-88cf-4fd5-ad50-8d76a1379a27
md"""
- For some special cases, the algorithm will return the hull, but in general some overestimation will occur.
"""

# ╔═╡ 0308378a-2d36-4c47-add7-0f53dc9011ca
md"""
## Preconditioning

- For the previous algorithms to work and be numerically stable, there are some requirements on the matrix $\mathbf{A}$.
- If those requirements are not met, one can try to precondition the problem with a *real* matrix $C$. that is apply the algorithms to the linear system

$
C\mathbf{Ax}=C\mathbf{b}$
- A particularly good choice is $C\approx A_c^{-1}$
"""

# ╔═╡ 9dba04bd-a8c3-44da-a1ab-855829c824ae
md"""
##
To understand the need for preconditioning, let us consider the following example.
"""

# ╔═╡ 8f68a5f0-4ec9-40d5-89be-cddf0c34a977
@bind N Slider(2:100; show_value=true, default=5)

# ╔═╡ 22abccab-5c93-4aa4-90e4-8209997a9a86
A1 = tril(fill(1..1, N, N))

# ╔═╡ e77a478d-0270-4120-a2d4-36efc455fb2e
b1 = [-2..2, fill(0..0, N-1)...]

# ╔═╡ 1ec87e5e-2f90-490b-a94f-05ea53d4cb1d
md"""
The correct solution is $[[-2, 2], [-2, 2], [0, 0], [0, 0], \ldots]$
"""

# ╔═╡ de8734a4-46b8-45ba-aff1-fa6a73e85d96
md"""
##
"""

# ╔═╡ e3bee176-9b9a-448f-b2d0-fa53cad35555
solve(A1, b1, GaussianElimination(), NoPrecondition())

# ╔═╡ 495e13e3-bf20-4541-ab6c-f6cb9572e9ed
solve(A1, b1, HansenBliekRohn(), NoPrecondition())

# ╔═╡ f223ac43-79cd-43be-b585-907a0bfb8f27
md"""
##
"""

# ╔═╡ acb54294-917d-41b0-b501-364b9228531b
solve(A1, b1, GaussianElimination(), InverseMidpoint())

# ╔═╡ ea231b9c-4cc9-4d45-b171-4c66536b02a1
solve(A1, b1, HansenBliekRohn(), InverseMidpoint())

# ╔═╡ 44813e6f-d7a2-47d6-9862-b9fd3e102ee2
md"""
##
### Downsides of preconditioning

- In general, the solution set of the preconditioned problem is **not** the solution set of of the original problem.
- Let us consider our previous 2D example.
"""

# ╔═╡ 1a1b8edf-bc1a-466e-b714-3eb71aa145ce
begin
	polytopes2 = solve(A, b, LinearOettliPrager(), InverseMidpoint())
	plot(UnionSetArray(polytopes2), ratio=1, label="preconditioned", legend=:right)
	plot!(UnionSetArray(polytopes), label="original", α=1) 
end

# ╔═╡ d8735a5e-8444-4724-9796-a59b37efd0e2
md"""
##
### Take-home lesson
- In general preconditioning is needed to achieve numerical stability.
- This may however enlarge the solution set.
- If preconditioning is not specified, the package performs some heuristic checks to decide a precondition strategy.
"""

# ╔═╡ bbd4086e-af3f-4743-b9fb-01559f567560
solve(A1, b1, GaussianElimination())

# ╔═╡ 1e9bf340-2dca-4602-a60c-b2ccdc9cfb85
md"""
## Conclusions

- [IntervalLinearAlgebra.jl](https://github.com/lucaferranti/IntervalLinearAlgebra.jl) offers (will offer) a toolbox to deal with interval linear systems.
- New features (possibly!) coming: 
  - Parametric interval systems
  - Eigenvalues computation
  - Determinant computation

- Ultimate goal: **Linear algebra done rigorously!**
- Download the slides: <https://github.com/lucaferranti/ILAjuliacon2021>
"""

# ╔═╡ Cell order:
# ╟─f403b050-e2e6-11eb-1c44-911de0d748bd
# ╟─0bc12fe2-4d08-4840-b651-b58df23f0277
# ╟─c9e05bd9-1d7b-41c0-99dc-29b1c1e20090
# ╟─c3a5fd72-9a9d-45f3-aaf3-bd61637d286c
# ╟─f2cff36e-9136-4eb1-acff-c97af31c2450
# ╟─286f40a2-840d-431b-a531-fe9e09dc55f1
# ╟─97e498f5-2e0c-4139-9774-86f58d59d511
# ╠═bb878728-d4ca-4c01-a4f7-bf1811a8f5c2
# ╠═a89ee6e9-565a-46d5-82db-bafb251a2200
# ╟─4848561e-2336-40fb-a5c6-8705da85c2ee
# ╟─5b111f71-4d41-43cf-a3cf-e2a6f4c2adba
# ╟─395b75e5-e76f-4f3a-b51f-57eaa3de2d4e
# ╟─f39f94f7-59ff-4dd0-93dd-4f9f158c673c
# ╟─3f5a6fed-6175-4459-8be2-9a330cedc938
# ╟─c223f95a-9885-48ae-b939-4f3fa54a7987
# ╟─32f0ff22-d42f-410f-a1a6-b8ea77b696ff
# ╟─a79e27c9-399c-44e1-ac79-f65480fe5178
# ╟─ba3a1ca6-04dc-45db-b686-a036519f2bc9
# ╟─38883c77-d506-4e81-9d25-0ac78323241a
# ╟─4549ca50-615c-4fba-afde-303ad231f648
# ╟─0f196850-f4b1-4830-9198-d72a281da05c
# ╟─d5e3856d-07e7-4e50-9e4d-c0afb7091d8b
# ╠═95edbbbf-2e35-446d-9eeb-a53f5c6b10ec
# ╟─4880c8a6-de84-4e55-be90-f73256b0888a
# ╟─75a66d2a-01b9-4d6a-8fed-8d8159ba6354
# ╟─cd558a80-2825-4953-870c-93551f60c10e
# ╟─c5aa282c-c82b-4521-994d-235cb7934c83
# ╟─32481e35-1a98-4b65-b925-264d51d24743
# ╟─58b8b1b8-93e7-46c0-b0a0-28532b152908
# ╟─faa03639-a431-4cbf-ac03-bbb33fb7b5b8
# ╠═5e0cf15b-21ab-45bb-81fa-ee962a46004f
# ╟─fc08d410-02fa-480b-a999-841476dfaaf6
# ╟─57a5fb2a-5cb2-412f-a1a0-5d27e0a5b8f7
# ╟─a70a6137-88cf-4fd5-ad50-8d76a1379a27
# ╟─0308378a-2d36-4c47-add7-0f53dc9011ca
# ╟─9dba04bd-a8c3-44da-a1ab-855829c824ae
# ╟─8f68a5f0-4ec9-40d5-89be-cddf0c34a977
# ╟─22abccab-5c93-4aa4-90e4-8209997a9a86
# ╟─e77a478d-0270-4120-a2d4-36efc455fb2e
# ╟─1ec87e5e-2f90-490b-a94f-05ea53d4cb1d
# ╟─de8734a4-46b8-45ba-aff1-fa6a73e85d96
# ╠═e3bee176-9b9a-448f-b2d0-fa53cad35555
# ╠═495e13e3-bf20-4541-ab6c-f6cb9572e9ed
# ╟─f223ac43-79cd-43be-b585-907a0bfb8f27
# ╠═acb54294-917d-41b0-b501-364b9228531b
# ╠═ea231b9c-4cc9-4d45-b171-4c66536b02a1
# ╟─44813e6f-d7a2-47d6-9862-b9fd3e102ee2
# ╟─1a1b8edf-bc1a-466e-b714-3eb71aa145ce
# ╟─d8735a5e-8444-4724-9796-a59b37efd0e2
# ╠═bbd4086e-af3f-4743-b9fb-01559f567560
# ╟─1e9bf340-2dca-4602-a60c-b2ccdc9cfb85
