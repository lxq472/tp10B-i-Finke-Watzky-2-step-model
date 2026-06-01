# Finke-Watzky 2-Step Model — ODE vs Analytical Solution

**Course:** Physics of Molecular Diseases — Week 1  
**Question 5:** Write code implementing the ODE described in the Finke-Watzky model. Test that your results agree with the analytical solution.

---

## Overview

The Finke-Watzky (F-W) 2-step model is a minimal "Ockham's razor" kinetic model for protein aggregation. Despite having only two parameters ($k_1$ and $k_2$), it fits 14 representative neurological protein aggregation datasets with $R^2 \geq 0.98$ — for amyloid- $\beta$, $\alpha$-synuclein, and polyglutamine disease proteins, across four different experimental methods and nine independent laboratories.

Its key advantage over the full Oosawa Master Equation approach is that it **deconvolutes nucleation from growth** in a simple, analytically tractable form.

---

## The Model

### Reactions (Scheme 2, Morris et al. 2008)

$$A \xrightarrow{k_1} B \qquad \text{slow continuous nucleation}$$

$$A + B \xrightarrow{k_2} 2B \qquad \text{fast autocatalytic surface growth}$$

where $A$ is the free monomer (unfolded/misfolded protein) and $B$ represents all aggregated species (oligomers, protofibrils, fibrils — lumped into one average species).

### ODE (eq. 1)

$$\frac{d[A]}{dt} = -k_1[A] - k_2[A][B]$$

with conservation law $[A] + [B] = [A]_0$, so $[B] = [A]_0 - [A]$.

The ODE is integrated numerically using `scipy.integrate.solve_ivp` (RK45).

### Analytical solution (eq. 2 & 3)

$$[A]_t = \frac{\dfrac{k_1}{k_2} + [A]_0}{1 + \dfrac{k_1}{k_2[A]_0} \exp\!\bigl((k_1 + k_2[A]_0)\,t\bigr)}$$

$$[B]_t = [A]_0 - [A]_t$$

This closed-form solution exists because the conservation law reduces the system to a single ODE in $[A]$, which is a Bernoulli equation with an exact solution.

---

## Parameters and their physical meaning

| Parameter | Units | Physical meaning | Observable |
|-----------|-------|-----------------|------------|
| $k_1$ | time$^{-1}$ | Nucleation rate constant | $k_1 \propto 1/t_\text{induction}$ — controls lag phase length |
| $k_2$ | conc$^{-1}$ time$^{-1}$ | Autocatalytic growth rate | $k_2[A]_0 \propto$ slope after lag — controls growth speed |
| $[A]_0$ | conc | Initial monomer concentration | Set experimentally |

The separation of $k_1$ and $k_2$ is the central result of the F-W model — something the earlier, more complete Oosawa-style models could not achieve in a simple, generally applicable way.

### Effect of parameters

| Change | Effect on $[B](t)$ |
|--------|-------------------|
| ↑ $k_1$ | Shorter lag phase ($t_\text{ind} \propto 1/k_1$); nucleation starts sooner |
| ↓ $k_1$ | Longer lag phase; delayed aggregation onset |
| ↑ $k_2$ | Steeper growth slope after lag; faster conversion to aggregates |
| ↓ $k_2$ | Slower growth; lag phase almost unchanged |

This is the key mechanistic insight: $k_1$ and $k_2$ control independent features of the sigmoidal curve and can therefore be independently extracted from experimental data.

---

## Validation of ODE vs Analytical

The numerical ODE solution (RK45, `rtol=1e-10`, `atol=1e-12`) agrees with the analytical solution to better than $2 \times 10^{-11}$ mM absolute error — a relative error of $~10^{-9}$ %. This confirms both implementations are correct and consistent.

The % aggregated curve ($[B](t)/[A]_0 \times 100$) reproduces Figure 1 of Morris et al. (2008) with the same representative parameters ($k_1 = 10^{-5}$ h$^{-1}$, $k_2 = 10^{-3}$ mM$^{-1}$h$^{-1}$).

---

## Comparison with the Oosawa / Master Equation approach

| Feature | Oosawa + Master Equation | Finke-Watzky |
|---------|--------------------------|--------------|
| Parameters | $k_n$, $k_a$, $k_d$, $k_f$, $n_c$ | $k_1$, $k_2$ only |
| Mechanistic detail | High — tracks f(j,t) | Low — B is a catch-all |
| Analytical solution | Only in limiting cases | Exact closed form |
| Fitting experimental data | Difficult | Straightforward |
| Deconvolutes nucleation from growth | Hard | Yes — directly |
| Predicts size distribution f(j) | Yes | No |

The F-W model is phenomenological — $k_1$ and $k_2$ are averages over potentially hundreds of elementary steps. This is its main limitation, but also what makes it broadly applicable.

---

## Limitations of the F-W model (from Morris et al. 2008)

1. **Oversimplified:** The real aggregation process involves hundreds to thousands of elementary steps; $k_1$ and $k_2$ are averages that hide important mechanistic detail.

2. **Size-independent rates:** The model assumes $k_1$ and $k_2$ are the same regardless of aggregate size, which is never exactly true and may hide size-dependent changes in growth rate.

3. **B is a catch-all:** All aggregate sizes (oligomers, protofibrils, mature fibrils) are lumped into B. Since smaller intermediates appear to be the more toxic species in neurological diseases, this is a significant weakness for understanding toxicity.

---

## Key References
- Morris, A. M., Watzky, M. A., Agar, J. N., & Finke, R. G. (2008). Fitting neurological protein aggregation kinetic data via a 2-step, minimal/"Ockham's Razor" model: the Finke-Watzky mechanism of nucleation followed by autocatalytic surface growth. *Biochemistry*, 47(8), 2413–2427.
- Prof. Ala Trusina, Lecture notes from the course "Physics of Molecular Diseases", Niels Bohr Institute, 2020
