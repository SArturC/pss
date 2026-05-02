# Robust Passive Systems via LMIs

This repository contains a database and tools for the analysis and synthesis of uncertain passive linear systems using Linear Matrix Inequalities (LMIs).

## 📦 Contents

- `data/`: Database of generated systems

## ⚙️ Requirements

- MATLAB / Octave
- YALMIP
- MOSEK (recommended)

## 🚀 Usage

Example:

```matlab
loaded = load('db_n2_m1_v3.mat');
f = fieldnames(loaded);
db = loaded.(f{1});

% Matrices
i   = 1 % i in [1, ..., 500]
Acl = db(i).Acl;
A   = db(i).A;
Bu  = db(i).Bu;
Bw  = db(i).Bw;
C   = db(i).C;
Du  = db(i).Du;
Dw  = db(i).Dw;
K   = db(i).K;
