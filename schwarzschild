(* Mathematica Notebook File *)
(* Created with Wolfram Language support in VS Code *)

(* ::Package:: *)

(* ::Title:: *)
(*Mathematica Notebook*)

(* GENERIC METRIC WORKFLOW - CHANGE PARAMETERS HERE *)

(* ==== CONFIGURATION SECTION ====*)
(* Define coordinates *)
coordinates = {t, r, theta, phi};

(* MODIFY THIS MATRIX FOR DIFFERENT METRICS *)
metricMatrix = {
  {-(1 - 2 M/r), 0, 0, 0},
  {0, (1 - 2 M/r)^(-1), 0, 0},
  {0, 0, r^2, 0},
  {0, 0, 0, r^2 Sin[theta]^2}
};

metricName = "Schwarzschild";

(* ==== AUTOMATIC CALCULATIONS ====*)
Print["===== ", metricName, " Metric Workflow =====\n"];
Print["Metric Tensor:"];
Print[Grid[metricMatrix, Frame -> All]];

(* Compute inverse metric *)
inverseMetricMatrix = Simplify[Inverse[metricMatrix]];
Print["\nInverse Metric Tensor:"];
Print[Grid[inverseMetricMatrix, Frame -> All]];

(* Christoffel symbols: Γ^λ_μν = (1/2) g^λσ (∂_μ g_σν + ∂_ν g_μσ - ∂_σ g_μν) *)
christoffel[lambda_, mu_, nu_] := Module[{sigma},
  (1/2) * Sum[
    inverseMetricMatrix[[lambda, sigma]] * (
      D[metricMatrix[[sigma, nu]], coordinates[[mu]]] +
      D[metricMatrix[[mu, sigma]], coordinates[[nu]]] -
      D[metricMatrix[[mu, nu]], coordinates[[sigma]]]
    ),
    {sigma, 1, 4}
  ]
];

Print["\nChristoffel Symbols (non-zero components):"];
Do[
  If[christoffel[i, j, k] != 0,
    Print["Γ^", i, "_{", j, k, "} = ", Simplify[christoffel[i, j, k]]]
  ],
  {i, 1, 4}, {j, 1, 4}, {k, 1, 4}
];

(* Riemann tensor: R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ *)
riemann[rho_, sigma_, mu_, nu_] := Module[{lambda},
  D[christoffel[rho, nu, sigma], coordinates[[mu]]] -
  D[christoffel[rho, mu, sigma], coordinates[[nu]]] +
  Sum[christoffel[rho, mu, lambda] * christoffel[lambda, nu, sigma], {lambda, 1, 4}] -
  Sum[christoffel[rho, nu, lambda] * christoffel[lambda, mu, sigma], {lambda, 1, 4}]
];

(* Ricci tensor: R_μν = R^λ_μλν *)
ricci[mu_, nu_] := Sum[riemann[lambda, mu, lambda, nu], {lambda, 1, 4}];

Print["\n\nRicci Tensor Components:"];
ricciMatrix = Table[Simplify[ricci[i, j]], {i, 1, 4}, {j, 1, 4}];
Print[Grid[ricciMatrix, Frame -> All]];

(* Ricci scalar: R = g^μν R_μν *)
ricciScalar = Sum[inverseMetricMatrix[[i, j]] * ricciMatrix[[i, j]], {i, 1, 4}, {j, 1, 4}];
Print["\nRicci Scalar:"];
Print[Simplify[ricciScalar]];

(* Einstein tensor: G_μν = R_μν - (1/2) g_μν R *)
Print["\nEinstein Tensor Components:"];
einsteinMatrix = Table[
  Simplify[ricciMatrix[[i, j]] - (1/2) * metricMatrix[[i, j]] * ricciScalar],
  {i, 1, 4}, {j, 1, 4}
];
Print[Grid[einsteinMatrix, Frame -> All]];

(* Kretschmann scalar: K = R_ρσμν R^ρσμν *)
Print["\nKretschmann Scalar (measure of curvature):"];
Print["Computation in progress..."];
