function [dY] = IDModelODE(t, Y, P)
% ODE for ID model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [IDH; QDFLocal; QDBLocal; IDF; IDB; QDF; QDB] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

ID = P.parameters.ID;
GC = P.parameters.GC;

%% Input
IDH       = Y(1);
QDFLocal  = Y(2);
QDBLocal  = Y(3);
IDF       = Y(4);
IDB       = Y(5);
QDF       = Y(6);
QDB       = Y(7); %NOTE: Re-removed the factor of 18...?

%% Variables
% Time dependent.
n = GetTimeIndex(t, P.results.tArray);
IBolus = P.results.IBolus(n); % [mU/min]

%% Computationl
dIDH = IBolus ...  % Injected subcutaneous insulin detemir.
    - ID.ka*IDH;   % Dissociation to mono/dimers in local interstitium.

dQDFLocal = ID.ka * IDH ...                    % Dissociation from hexamers.
    + [ID.kd2*QDBLocal - ID.kd1*QDFLocal] ...  % Binding exchange.
    - (ID.kb+ID.kdi) * QDFLocal;               % Transport to plasma / degradation.

dQDBLocal = ID.kd1 * QDFLocal...  % Binding exchange.
    - ID.kd2 * QDBLocal;

dIDF = ID.kb/GC.VI * QDFLocal ...    % Transport from local interstitium.
    + [ID.kd2*IDB - ID.kd1*IDF] ...  % Binding exchange.
    - ID.nDI*(IDF - QDF) ...         % Exchange with interstitium.
    - (ID.nDL + ID.nDK) * IDF;       % Hepatic / renal clearance.

dIDB = ID.kd1*IDF ...  % Binding exchange.
    - ID.kd2*IDB;

dQDF =  ID.nDI*(IDF - QDF) ...      % Exchange with plasma.
    + [ID.kd2*QDB - ID.kd1*QDF]...  % Binding exchange.
    - ID.nDC*QDF;                   % Degradation.

dQDB = ID.kd1*QDF ...  % Binding exchange.
    - ID.kd2*QDB;

%% Output
dY = [dIDH;
    dQDFLocal;
    dQDBLocal;
    dIDF;
    dIDB;
    dQDF;
    dQDB];

end
