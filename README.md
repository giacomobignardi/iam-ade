# Assortative mating can mask non-additive genetic contributions to individual differences

This repository contains the full implementation of the **iAM-ADE model**. It includes:

- Code to execute the iAM-ADE model
- Scripts to reproduce the simulation study
- Code and resources for the analyses presented in:

> Bignardi, G., Sunde, H. F., Bruins, S., Balbona, J. V., Fisher, S. E., & Boomsma, D. I. (*in revision*). Assortative mating can mask non-additive genetic contributions to individual differences.

---

## Repository Structure

```
.
├── R/                          # Scripts for the analysis
│   └── functions/
│       ├── acde.fun            # Classic twin model
│       ├── iam.date.fun        # Extended twin family model
│       ├── iam.date.fun.helper # iAM-ADE model specification
│       ├── iam.date.fun.plot   # Plotting utilities for iAM-ADE outputs
│       └── sim.estd.fun        # Simulate twin, sibling, and partner data
│                               # under indirect assortative mating
│                               # (and other mechanisms of family resemblance)
├── SI/                         # Supporting information
├── LICENSE
├── README.md
└── .gitignore
```

---

## License

Copyright (c) 2026 [Authors]. This code is made available for review purposes only. No use, modification, or distribution is permitted without explicit written consent from the authors. A permissive open-source license will be issued upon publication.
