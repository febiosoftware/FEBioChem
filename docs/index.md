---
date:
  created: 2026-1-18
---

![FEBioChem logo](./figs/FEBioChem_logo_high_res.png){ style="display:block; margin:0 auto; width:50%;" }

## Overview

**FEBioChem** is an FEBio plugin for modeling **reaction–diffusion** and **reaction–diffusion–convection** processes in non-deformable mixtures. It is designed for systems involving multiple interacting chemical species, where transport occurs through diffusion and advection, and species may also participate in complex reaction networks. These species can exist either as **mobile solutes** carried by the solvent or as **solid-bound constituents** attached to the mixture matrix. By assuming the mixture is non-deformable, FEBioChem eliminates the need to solve for solid mechanics, leading to **simplified formulations**, **reduced computational cost**, and **improved performance**—while still capturing the essential physics of chemically reactive transport.

---

## Getting Started

- 📘 **[Tutorial](febiochem_tutorial1.md)**  
  Step-by-step guide using FEBioChem with FEBio Studio  

- 📖 **[User Manual](febiochem_user_manual.md)**  
  Detailed reference for input file configuration  

- 🧠 **[Theory](febiochem_theory.md)**  
  Background on governing equations and implementation details 

## Resources
The easiest way to get to the FEBioChem plugin is to use the Plugin Repo in FEBio Studio. Just open FEBio Studio, use the menu **FEBio → Plugin Repository** to open the plugin repo, then find and download the FEBio Chem plugin. 

Some additional resources are linked below. 

:material-github: [FEBioChem](https://github.com/febiosoftware/FEBioChem) on Github.

:material-link-variant: [FEBioChem](https://repo.febio.org/pluginRepo/Plugins?id=2){ target="_blank" rel="noopener noreferrer" } on the Plugin Repo website.

:material-link-variant: The [FEBio](https://febio.org/){ target="_blank" rel="noopener noreferrer" } website 