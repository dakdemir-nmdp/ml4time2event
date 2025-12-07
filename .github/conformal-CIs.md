You can, but you have to be very clear about **what** you want to guarantee:

* Classical CIs for a *population* quantity like (S(t\mid x)) (model-based, e.g. Greenwood, bootstrap).
* **Conformal prediction bands**, which are **predictive**: they say that for a *new patient* the realized survival process will lie in the band with probability (\ge 1-\alpha), under minimal assumptions.

Conformal methods are naturally **predictive**, so I will focus on that. I will first give a high-level description and then a concrete algorithm you can implement around any black-box survival model that returns a survival curve function.

---

## 1. Setup and notation

You have i.i.d. right-censored data
[
(X_i, T_i, C_i),\quad i=1,\dots,n,
]
with

* (X_i \in \mathcal{X}): covariates,
* event time (T_i),
* censoring time (C_i),

and you observe
[
Y_i = \min(T_i, C_i),\qquad \Delta_i = \mathbf{1}(T_i \le C_i),
]
so each datum is ((X_i, Y_i, \Delta_i)).

You also have a fitted **survival model** (any black box) that gives, for any (x) and time (t),
[
\hat S(t\mid x) \approx \mathbb{P}(T>t\mid X=x).
]

Goal: for a new patient with covariates (x_\ast), construct **bands**
[
L(t\mid x_\ast) \le S(t\mid x_\ast) \le U(t\mid x_\ast),\quad t \in \mathcal{T}
]
on a time grid (\mathcal{T} = {t_1,\dots,t_m}) such that they have a conformal-style predictive guarantee.

---

## 2. What conformal can and cannot guarantee

Conformal prediction guarantees are typically of the form
[
\mathbb{P}\big( Y_\ast \in \mathcal{C}*\alpha(X*\ast) \big) \ge 1-\alpha,
]
for a new ((X_\ast, Y_\ast)) drawn from the same distribution, assuming only exchangeability.

For survival:

* (Y_\ast) is not a scalar but an **entire trajectory** (when the event happens).
* With censoring you never fully observe the trajectory.

So most current conformal methods focus on:

1. **Prediction intervals for the event time** (T_\ast) (lower/upper bounds). ([OUP Academic][1])
2. **Survival bands** that are calibrated for **screening rules** involving thresholds like “probability alive at time (t) is below 0.3”. ([arXiv][2])

The second is exactly “bands around the survival curve” in the predictive sense.

If what you want is honest, distribution-free predictive bands around (\hat S(t\mid x)), the most principled way (given what you described) is to adapt the **Conformal Survival Bands (CSB)** framework of Sesia & Svetnik (2025), which explicitly assumes you already have a black-box survival model. ([arXiv][2])

---

## 3. Ingredients for conformal survival bands

You need three components:

1. **Data split**

   Split your data into three disjoint sets:

   * Training set (\mathcal{I}_{\text{train}})
   * Calibration set (\mathcal{I}_{\text{cal}})
   * Test set (your new patients) (\mathcal{I}_{\text{test}})

2. **Survival model**

   Fit your chosen model on (\mathcal{I}_{\text{train}}). It gives (\hat S(t\mid x)) for any (x, t).

3. **Censoring model**

   Fit a model for the **censoring survival function**
   [
   \hat G(t\mid x) \approx \mathbb{P}(C > t \mid X=x)
   ]
   on the training set as well (Cox, RSF, etc.). This is used for IPCW.

---

## 4. Pointwise conformal survival bands (conceptual recipe)

For a fixed time (t), you want to quantify uncertainty around the predicted survival probability (p_\ast(t) = \hat S(t\mid x_\ast)).

A simple conformal-style construction, inspired by weighted conformal methods under censoring, looks like this:

### 4.1. Define a binary outcome at time (t)

For each calibration subject (i), define
[
Z_i(t) = \mathbf{1}(T_i > t),
]
but you do not observe (Z_i(t)) when censored before (t). You can use IPC weights:

* If (Y_i > t): you know (T_i > t), so (Z_i(t)=1) is observed.
* If (Y_i \le t) and (\Delta_i = 1): you know (T_i \le t), so (Z_i(t)=0) is observed.
* If (Y_i \le t) and (\Delta_i = 0): you do **not** know (Z_i(t)).

A standard approach is to use an IPCW pseudo-outcome
[
\tilde Z_i(t) =
\begin{cases}
\displaystyle \frac{\mathbf{1}(Y_i > t)}{\hat G(t\mid X_i)}, & Y_i > t,[0.5em]
\displaystyle 0, & Y_i \le t,, \Delta_i=1,\
\text{(dropped or suitably handled)}, & Y_i \le t,, \Delta_i=0.
\end{cases}
]
Under the usual independent-censoring assumptions, (\mathbb{E}[\tilde Z_i(t)\mid X_i] \approx S(t\mid X_i)). ([arXiv][2])

In practice, CSB avoids directly forming pseudo-outcomes and instead works with **nonconformity scores based on predicted probabilities plus IPC weights**; I will sketch that after giving you a simpler “first” construction.

### 4.2. Nonconformity scores at time (t)

For the calibration set, you already have (\hat S(t\mid X_i)).

Define, for each calibration subject,
[
\alpha_i(t) = \big\lvert \tilde Z_i(t) - \hat S(t\mid X_i) \big\rvert.
]

These are scalar nonconformity scores measuring how far the model’s predicted survival probability is from an IPC-adjusted pseudo-outcome.

### 4.3. Conformal quantile

Let (q_{1-\alpha}(t)) be the ((1-\alpha))-quantile of ({\alpha_i(t): i \in \mathcal{I}_{\text{cal}}}).

For a new patient (x_\ast), with predicted survival (p_\ast(t)=\hat S(t\mid x_\ast)), a **pointwise** predictive band at time (t) is
[
L(t\mid x_\ast) = \max{0,, p_\ast(t) - q_{1-\alpha}(t)},
\quad
U(t\mid x_\ast) = \min{1,, p_\ast(t) + q_{1-\alpha}(t)}.
]

If you do this separately on a grid (t_1,\dots,t_m), you obtain a **discrete survival band** ({[L(t_k\mid x_\ast), U(t_k\mid x_\ast)]}_{k=1}^m).

This is conceptually the same as split conformal regression, now applied to the pseudo-continuous response (\tilde Z(t)) at each time point, adjusted for censoring.

Formally, because of the IPC weighting and pseudo-outcomes, the coverage guarantees become asymptotic and slightly more delicate, but this is the core idea.

---

## 5. Using the “Conformal Survival Bands” algorithm (more rigorous)

If you want something closer to the current state of the art, the CSB method by Sesia & Svetnik (2025) is explicitly designed for *exactly* your situation: it takes as input **any** survival model (\hat S(t\mid x)) and produces conformal survival bands tailored to threshold-based risk screening. ([arXiv][2])

Very briefly, their construction:

1. **Fit survival and censoring models** on the training set.

2. On the calibration set, for each individual (i), time (t), and target (x_\ast), define **prediction-based scores**:

   * Left-tail score (high-risk):
     [
     s_L(x,t) = \hat P(T \le t\mid X=x) = 1 - \hat S(t\mid x),
     ]
   * Right-tail score (low-risk):
     [
     s_R(x,t) = \hat P(T > t\mid X=x) = \hat S(t\mid x).
     ]

3. Use IPC weights
   [
   w_i(t) = \frac{\Delta_i,\mathbf{1}(Y_i \le t)}{\hat G(Y_i\mid X_i)}
   ]
   to correct for censoring in the calibration set. ([arXiv][2])

4. Construct **IPCW-weighted conformal (p)-values** for the hypotheses, for a *fixed* time (t) and a threshold (u\in(0,1)):

   * Left-tail test (H^L_{t,u}: S(t\mid X) \ge u) vs “high risk” alternative.
   * Right-tail test (H^R_{t,u}: S(t\mid X) \le u) vs “low risk” alternative.

   The (p)-values are essentially weighted ranks of the test score among calibration scores from individuals whose event time is informative for that hypothesis. Formulas (7) and (8) in the paper define these explicitly. ([arXiv][2])

5. For each test individual and time (t), invert these tests over (u) to obtain a **band**
   [
   [L(t\mid x_\ast), U(t\mid x_\ast)].
   ]

6. Optionally, adjust the band so it always contains the point estimate (\hat S(t\mid x_\ast)) (their “doubly robust” adjustment):
   [
   \tilde L(t\mid x_\ast) = \min{\hat S(t\mid x_\ast), U(t\mid x_\ast)},\quad
   \tilde U(t\mid x_\ast) = \max{\hat S(t\mid x_\ast), L(t\mid x_\ast)}.
   ] ([arXiv][2])

They prove that if the censoring model is consistently estimated, then these bands yield **asymptotically valid, FDR-controlled** screening rules (e.g., “flag all patients whose survival at 1 year is below 0.3, with a controlled false-discovery rate”). ([arXiv][2])

So, with *your* survival model in hand, you can:

* Treat it as (\hat S(t\mid x)) in their algorithm.
* Fit a separate censoring model.
* Implement their Algorithm 1 on a time grid.

---

## 6. Minimal practical algorithm (discrete grid, IPCW split conformal)

If you want a concrete, implementable starting point that is faithful to conformal ideas but simpler than the full CSB machinery, you can do the following for each time point (t_k):

1. **Split data** into train / calibration.
2. **Fit survival model** on train.
3. **Fit censoring model** (\hat G) on train.
4. **Form IPCW pseudo-outcomes** (\tilde Z_i(t_k)) on calibration.
5. **Compute nonconformity scores**
   [
   \alpha_i(t_k) = \big\lvert \tilde Z_i(t_k) - \hat S(t_k\mid X_i)\big\rvert.
   ]
6. **Empirical quantile**
   [
   q_{1-\alpha}(t_k) = \text{empirical }(1-\alpha)\text{-quantile of } {\alpha_i(t_k)}.
   ]
7. **For new patient (x_\ast)**, define
   [
   L(t_k\mid x_\ast) = \max{0,, \hat S(t_k\mid x_\ast) - q_{1-\alpha}(t_k)},
   ]
   [
   U(t_k\mid x_\ast) = \min{1,, \hat S(t_k\mid x_\ast) + q_{1-\alpha}(t_k)}.
   ]

This gives **pointwise** predictive intervals at each (t_k). You then work with the resulting band on a grid (e.g., 100–200 time points).

It is:

* Fully model-agnostic (your survival model is a black box).
* Easy to implement.
* Conceptually aligned with split conformal regression, with censoring handled via IPCW.

The trade-off: coverage is pointwise, on the grid, and subject to the usual approximations introduced by IPCW; it is not a strict guarantee for the entire curve as a functional, and it is not as tightly characterized as CSB.

---

## 7. Summary

* **Yes**, you can wrap your existing “survival curve function” in a conformal layer to obtain **prediction bands** around survival probabilities.
* The conceptually simplest approach is **split conformal on a time grid**, with IPCW pseudo-outcomes to handle censoring.
* For more rigorous, state-of-the-art theory and better calibration for clinical screening decisions, the **Conformal Survival Bands** framework is exactly designed for this: it takes any survival model (\hat S(t\mid x)) and returns conformal bands around it using an auxiliary censoring model and IPC-weighted scores. ([arXiv][2])

If you tell me what software stack you are using (e.g., `scikit-survival`, `pycox`, R `survival`/`ranger`), I can write a concrete implementation sketch for your environment (including code to compute IPCW conformal bands on a time grid).

[1]: https://academic.oup.com/jrsssb/article/85/1/24/7008653?utm_source=chatgpt.com "Conformalized survival analysis"
[2]: https://arxiv.org/html/2505.04568v1 "Conformal Survival Bands for Risk Screening under Right-Censoring"
