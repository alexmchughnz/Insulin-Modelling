function P = TagPatientCode(P, tag)

    P.patientCode = strjoin([P.patientCode, "("+tag+")"]);
end