# influential-groups

Code associated with the paper:<br>
**Influential groups for seeding and sustaining nonlinear contagion in heterogeneous hypergraphs**.

# Reference

If you use this code, please consider citing:

"[_Influential groups for seeding and sustaining nonlinear contagion in heterogeneous hypergraphs_](https://www.nature.com/articles/s42005-021-00788-w)" <br>
Guillaume St-Onge, Iacopo Iacopini, Vito Latora, Alain Barrat, Giovanni Petri, Antoine Allard, Laurent Hébert-Dufresne <br>
Commun. Phys., **5**, 25 (2022)


# Requirements

The scripts make use of 3 different submodules:
- [gcm](https://github.com/gstonge/gcm) for approximate master equations
- [horgg](https://github.com/gstonge/horgg) for the synthetic generation of hypergraphs
- [schon](https://github.com/gstonge/schon) for the simulation of nonlinear contagions on hypergraphs
Please refer to each submodule for their own requirements.

We also make use of standard python packages `numpy`, `matplotlib`,
and many others.

# Installation

First clone the repository and the submodules using

```
git clone --recurse-submodules git@github.com:gstonge/influential-groups.git
```

or

```
git clone --recurse-submodules https://github.com/gstonge/influential-groups.git
```

Depending on the scripts you want to use, you might need to install `gcm`, `horgg`, and/or `schon`.
Please refer to each submodule installation instructions.


