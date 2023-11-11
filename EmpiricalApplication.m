%*********************************************************************
%    In this program, we conduct the CUSUM and QS tests by 
%    Zhang, Wu and Wu (Statistical Research ) to log return data of
%    exchange rate between dollar and rouble during the Russia-Ukraine War.
%
%    by Ruike Wu @ Xiamen University
% *******************************************************************

clear;clc;
load('Dollar_Rouble_log_return.mat');
y = 100.*log_return1';
CUSUM_setval = [];
QS_setval = [];
% Conduct our test using AR(1) Plug-in procedure
[CUSUM_setval(1) ,QS_setval(1)] = Vol_structural_change_test(y,-1);

