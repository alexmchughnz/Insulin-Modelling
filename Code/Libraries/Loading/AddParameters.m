function P = AddParameters(P)
P.parameters.CP = CPParameters(P);
P.parameters.GC = GCParameters(P);
P.parameters.GI = GIParameters(P);
if P.data.IType == "detemir"
    P.parameters.ID = IDParameters(P);
elseif P.data.IDelivery == "subcutaneous"
    P.parameters.SC = SCParameters(P);
end
end
