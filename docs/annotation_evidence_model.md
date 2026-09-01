# Annotation Evidence Model

The model follows the supplied Shin-MassBank proposal and keeps three dimensions
separate:

1. Annotation level: `L1`, `L2`, `L3`, `L4`, `L5`.
2. Claim specificity for Level 3: `SP`, `CP`, `MO`, `CL`.
3. Evidence inventory: `RS`, `FM`, `SL`, `DF`, `RT`, `IM`, `IS`, `MN`, `CX`.

Examples include `L1[RS,RT,SL,FM]`, `L2[SL,RT,FM]`, and
`L3-CP[FM,DF,MN,CX]`. A claim is not promoted solely because it has many evidence
tags. A versioned `criteria_set` defines which evidence and thresholds are needed
for a study, and measured values such as mass error, RT error, dot products, CCS
error, edge score, or q-like values are retained with each assertion.

