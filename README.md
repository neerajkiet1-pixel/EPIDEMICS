# EPIDEMICS
Epidemics spread over n weeks
# Epidemic Spread Over N Weeks: Project Briefing

## 1. What the problem means in simple terms

Picture a population of 1,000 people. Each week, every person is in one of three states:

- **S (Susceptible):** healthy, can catch the disease
- **I (Infected):** currently sick
- **R (Recovered):** immune, at least for now

Each week, some people move between states, and the probabilities of those moves don't change. We collect those probabilities in a 3×3 matrix **A**. If you know today's split (say 99% S, 1% I, 0% R), multiplying by A gives next week's split. Multiplying 52 times gives the split after one year.

So the project asks three things:
1. How do we get week 52 without doing 52 separate multiplications? (Cayley-Hamilton)
2. What does the whole year look like, week by week?
3. Where does the epidemic settle in the long run? (steady state)

**One honest caveat to know now:** real SIR models are nonlinear, because infection rate depends on how many people are already infected. A fixed matrix A assumes constant probabilities, so this is a *linear Markov-chain approximation*. It's valid for the early phase, for a fixed intervention level, or as a discrete-time teaching model. Your report should say this. Examiners like it when you show you know the model's limits.

---

## 2. The exact mathematics

### 2.1 Cayley-Hamilton theorem
Every square matrix satisfies its own characteristic equation. For a 3×3 matrix:

p(λ) = det(λI − A) = λ³ − c₂λ² + c₁λ − c₀ = 0

where
- c₂ = trace(A)
- c₁ = sum of the principal 2×2 minors
- c₀ = det(A)

Cayley-Hamilton says that substituting A for λ gives

**A³ = c₂A² − c₁A + c₀I**

So A³ is a combination of A², A and I. By extension, **every higher power Aⁿ can be written using only I, A and A²**:

**Aⁿ = a₂(n)·A² + a₁(n)·A + a₀(n)·I**

The whole problem reduces to finding three scalar coefficients for n = 52.

### 2.2 Finding the coefficients
Divide the polynomial λⁿ by p(λ):

λⁿ = q(λ)·p(λ) + r(λ), where r(λ) = a₂λ² + a₁λ + a₀

Since p(A) = 0, the quotient term vanishes and Aⁿ = r(A).

There are two ways to find r(λ):

- **Eigenvalue method (distinct roots):** For each eigenvalue λᵢ, λᵢⁿ = r(λᵢ). That gives 3 linear equations in a₀, a₁, a₂ (a Vandermonde system).
- **Polynomial reduction (always works, even with repeated roots):** Compute λⁿ mod p(λ) by repeated squaring, reducing the degree each step.

### 2.3 Why it's efficient
Computing A⁵² naively takes 51 matrix multiplications. With Cayley-Hamilton, you work with *scalar polynomials of degree ≤ 2* and form one final combination. For a 3×3 matrix the savings are modest in absolute terms, but the idea scales: for an n×n matrix, you never need powers above Aⁿ⁻¹. Be upfront about this in your presentation, and verify your result against NumPy's `matrix_power`.

---

## 3. Input variables

| Symbol | Meaning | Typical range |
|---|---|---|
| β | Weekly probability a susceptible person gets infected | 0.05 to 0.5 |
| γ | Weekly probability an infected person recovers | 0.2 to 0.7 (≈ 1/duration in weeks) |
| δ | Weekly probability a recovered person loses immunity (back to S) | 0 to 0.1 |
| x₀ = (S₀, I₀, R₀) | Initial population fractions, summing to 1 | e.g. (0.99, 0.01, 0) |
| n | Number of weeks | 52 (user adjustable) |

---

## 4. Forming the model

Use the state vector **xₖ = [Sₖ, Iₖ, Rₖ]ᵀ** and the rule **xₖ₊₁ = A·xₖ**, so **xₙ = Aⁿ·x₀**.

I'd recommend the SIRS version (with waning immunity). Columns are "from" states, rows are "to" states:

```
        from S    from I    from R
to S  [ 1−β        0         δ    ]
to I  [  β        1−γ        0    ]
to R  [  0         γ        1−δ   ]
```

Every **column sums to 1** (column-stochastic), because everyone leaves a state to go somewhere, including staying put.

**Why not plain SIR (δ = 0)?** Then A is lower triangular with eigenvalues 1−β, 1−γ and 1. The only steady state is (0, 0, 1), where everyone ends up recovered. That's correct but dull. Making δ an optional parameter lets the demo show both cases.

**Characteristic polynomial pieces for the SIRS matrix:**
- trace = 3 − β − γ − δ
- c₁ = (1−β)(1−γ) + (1−β)(1−δ) + (1−γ)(1−δ)
- det = (1−β)(1−γ)(1−δ) + βγδ

Because A is column-stochastic, **λ = 1 is always an eigenvalue**. That fact drives the steady state.

**Steady state:** Solve A·s = s with S + I + R = 1. For δ > 0 this gives

s : i : r = 1/β : 1/γ : 1/δ

So S* = (1/β) / (1/β + 1/γ + 1/δ), and similarly for I* and R*. The other two eigenvalues have magnitude < 1 (under normal parameter choices), so their terms die out as n grows and the system converges to this vector regardless of x₀. This is a good point to explain in your report: **Aⁿ → a rank-1 matrix whose columns are all s**.

---

## 5. Python algorithm (to implement later)

1. **Input:** β, γ, δ, x₀, n. Validate that probabilities are in [0,1] and columns sum to 1.
2. **Build A.**
3. **Characteristic polynomial:** compute trace, c₁ and det, giving p(λ).
4. **Reduce λⁿ mod p(λ):** use repeated squaring on coefficient triples (52 = 110100 in binary, so about 6 squarings plus a few multiplies), reducing degree 3 → 2 each time using p(λ). This gives (a₂, a₁, a₀).
5. **Assemble Aⁿ = a₂A² + a₁A + a₀I.**
6. **Evolution:** compute xₖ for k = 0..52. You can use the coefficient recurrence, or call the Cayley-Hamilton routine for each k.
7. **Verify:** compare with `numpy.linalg.matrix_power(A, 52)` and report the maximum error.
8. **Steady state:** solve (A − I)s = 0 with normalization, and cross-check using the eigenvector for λ = 1 and the limit of Aⁿ for large n.
9. **Analysis extras:** eigenvalues, convergence rate (|λ₂|), peak infection week and peak value.

---

## 6. Final demo application

- **Inputs:** sliders or fields for β, γ, δ, initial infected fraction, and number of weeks.
- **Outputs:**
  - the matrix A, its characteristic polynomial and eigenvalues
  - the computed A⁵², shown next to NumPy's result with the error
  - a line chart of S, I and R over weeks 0 to 52
  - a table of weekly values
  - the steady-state vector, marked as a dashed line on the chart
  - key stats: peak infection week, week when the system is within 1% of steady state
- **Comparison mode:** timing of Cayley-Hamilton vs. naive repeated multiplication vs. NumPy.
- **Scenario presets:** "No immunity loss" (ends all-recovered), "Waning immunity" (settles into an endemic level), and "Lockdown" (lower β).
- **Suggested format:** a Streamlit app or a Jupyter notebook with interactive widgets.

---

**Next step:** if you're comfortable with this, I can walk through a fully worked numerical example by hand (e.g. β = 0.3, γ = 0.4, δ = 0.05) so you can check your code against it, then we move to the code.

You also wrote "Read me file." Do you want this explanation turned into a `README.md` for your project folder?
