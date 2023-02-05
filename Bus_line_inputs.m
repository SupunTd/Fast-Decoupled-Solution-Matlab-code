% BUS TYPES: TYPE-1=SLACK, TYPE-2=PV, TYPE-3=PQ

function [linedata,busdata,tol,acceleration_factor,npv,nlb,lpv] = Bus_line_inputs()     
%         |  From |  To   |   R     |   X     |   B/2     |  a   |
%         |  Bus  | Bus   |         |         |           |      |
linedata = [ 1      2       0.00    0.038     0.00000   1.0;
             1      3       0.000    0.031     0.0000 1.0;
             2      3      0.000    0.028     0.00000   1.0;
             ];

%         |Bus | Type | Vsp | theta |  PGi   |  QGi   |  PLi   |   QLi  | Qmin | Qmax | B_shunt  |
busdata = [ 1     1    1.04     0      0.0      0.0      0.0     0.000     0     0.0    0.0;
             2    2    1.01     0      3.0     0.0      0.0     0.0    0.00   +0.00      0.0;
            3     3    1.00     0      0.0      0.0      4.0     1.5     0     0.0    0.0;];
         
nlb=1;    
npv=1;        
lpv=2; 
tol=input('Enter the tolerance vALUE :- ');        
acceleration_factor=0;
