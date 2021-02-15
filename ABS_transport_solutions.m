
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%  ABC_TRANSPORT_SOLUTIONS.M                                               %
%                                                                          %
% Solve for analytical solution of uptake rate for ABC transport model.    %
%                                                                          %
%                                                                          %
%       Noele Norris                                                       %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


syms k0f k0r k1f k1r k2 k3
syms BP SBP TSBP TBP T S BP_total T_total


eqn1 = TSBP == (k3/(k2+k3))*T_total*SBP/(SBP + k3*(k2)/(k2+k3)/k1f);
eqn2 = SBP == S*BP/(k0r/k0f + k1f*T/k0f);
eqn3 = BP_total == BP + SBP + (1+k2/k3)*TSBP;
eqn4 = T_total == T + (1+k2/k3)*TSBP;

solutions = solve([eqn1, eqn2, eqn3, eqn4], {BP SBP TSBP TBP T})

TSBP = simplify(solutions.TSBP)
SBP = simplify(solutions.SBP)


