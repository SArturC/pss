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
load('data/database_unificada.mat');

system = db(1);
