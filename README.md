# Logic cloning based approximate circuits 

## Table of Contents
- [Project Overview](#project-overview)
- [Folder Structure](#folder-structure)
  - [1. Multiplication](#multiplication)
  - [2. RTL](#rtl)
  - [3. MIMO](#mimo)


## Project Overview
Please refer to following publication:
>Abhinav Kulkarni, Messaoud Ahmed Ouameur, Daniel Massicotte, Logic cloning based approximate signed multiplication circuits for FPGA,
Microelectronics Journal,2024, 106135,ISSN 0026-2692,(Online: https://www.sciencedirect.com/science/article/pii/S0026269224000478)

The project provides code for the work in above publication. Description of the project is provided as follows. 

## Folder Structure

### 1. Multiplication

The Multiplication folder contains C code implementations for various multiplication circuits. The `app.c` file holds the code, and `app.h` includes necessary definitions. This project is configurable in any Integrated Development Environment (IDE) like Codeblocks with the GCC compiler. It allows the configuration of multiplication circuits for different bit widths. Configure the IDE to generate `.so` static library files from `app.c`. These static library files are then utilized in Python to simulate multiplication operations with approximate multiplication circuits. The execution of this project is essential in a Linux environment.

### 2. RTL

The RTL folder includes VHDL code for the RTL implementation of signed multiplication. Separate files are provided for accurate and LC (Logic Cloning) variants. The accurate multiplication utilizes Baugh Wooley (BW) and Radix-2 Booth schemes, with code in files `accurate_bw.vhd` and `accurate_booth.vhd`. LC methodology is employed for approximate multiplication, and the corresponding code is found in files `LC_bw.vhd` and `LC_Booth.vhd`. These multiplication implementations are instantiated in the `top.vhd` file as the top component. The `top_multiple.vhd` file contains multiple instances of the top component, serving as the apex file in the file hierarchy. Utilize these files to initiate an RTL project.

### 3. MIMO 

The MIMO folder focuses on signal fidelity evaluation, measured in terms of Symbol Error Rate (SER), for MIMO uplink detection using Python. The `MIMO_SER.py` file assesses the SER of MIMO uplink detection employing Zero Forcing (ZF) and Minimum Mean Squared Error (MMSE) detection techniques. These techniques are implemented with approximate multiplication circuits. The `floatfix.py` file handles `c_types` operations, incorporating approximate multiplication and matrix inversion using the Gauss-Jordan method, leveraging approximate multiplication circuits. Note that the execution of this code requires a Linux environment.
