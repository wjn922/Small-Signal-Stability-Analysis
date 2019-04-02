# Small-Signal-Stability-Analysis
* written in C++
* This project refers to the SSS of a large-scale AC-DC power system.

## How to use
open the project SSSA_V2.sln, F5 run

## Characters
#### 1、Currently, the project includes the generator model of card MF and MG， and one excitor model of card EA.
#### 2、Input file need: 
* (1) Node table
* (2) Admittance matrix
* (3) Power flow
* (4) BPA .swi card
#### 3、Output file
* The linear model in the form of matrices (J,T)

## Files
#### 1、SSSA_V2 outputs the matrices (J,T) in the form of two-dimension array
#### 2、SSSA_V2 - sparse uses the sparse technique and outputs the matrices (J,T) in the form of CRS
