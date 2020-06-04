function [R] = IDModel(R, O, V)
% Function for ID model forward simulation.
% INPUTS:
%   R  - results struct, must have tArray
%   O  - ode45 options
%   V  - parameter variants
% OUTPUT:
%   R  - results struct updated with model results 

global ID

% Set up initial conditions.
Y0 = [ID.ISC0;
        ID.QDFLocal0;
        ID.QDBLocal0;
        ID.IDF0;
        ID.IDB0;
        ID.QDF0;
        ID.QDB0];
  
% Forward simulate.
[~, Y] = ode45(@IDModelODE, R.tArray, Y0, O);  

% Store results.
R.IDH      = Y(:,1);
R.QDFLocal = Y(:,2);
R.QDBLocal = Y(:,3);
R.IDF      = Y(:,4);
R.IDB      = Y(:,5);
R.QDF      = Y(:,6);
R.QDB      = Y(:,7);

end


function [dY] = IDModelODE(t, Y)
% ODE for ID model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [IDH; QDFLocal; QDBLocal; IDF; IDB; QDF; QDB] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global ID

%% Input
IDH       = Y(1);
QDFLocal  = Y(2);
QDBLocal  = Y(3);
IDF       = Y(4);
IDB       = Y(5);
QDF       = Y(6);
QDB       = Y(7);

%% Variables
% Time dependent.
n = (1 + floor(t));
IBolus = P.IBolus(n);

%% Computation
dIDH = -ID.ka*IDH + IBolus;
dQDFLocal = ID.ka*IDH - QDFLocal*(ID.kb + ID.kdi) ...
                - (ID.kd1*QDFLocal - ID.kd2*QDBLocal);
dQDBLocal = ID.kd1*QDFLocal - ID.kd2*QDBLocal;
dIDF = ID.kb/GC.VI(P)*QDFLocal - IDF*(ID.nDL + ID.nK) ...
           - ID.nDI*(IDF - QDF) - (ID.kd1*QDF - ID.kd2*QDB);
dIDB = ID.kd1*IDF - ID.kd2*IDB;
dQDF = -ID.nDC*QDF + ID.nDI*(IDF - QDF) - (ID.kd1*QDF - ID.kd2*QDB);
dQDB = ID.kd1*QDF - ID.kd2*QDB;

%% Output
dY = [dIDH;
      dQDFLocal;
      dQDBLocal;
      dIDF;
      dIDB;
      dQDF;
      dQDB];

end
