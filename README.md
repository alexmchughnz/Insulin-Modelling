A basic tutorial of glycaemic modelling or the STAR model.

For each patient, the code makes a patient struct "P" and runs it through `SimpleSim.m`, where most of the work is done.

Trial data is entered in `MakeTemplate.m`. Either load from spreadsheet (like the example `MakeOGTTLui.m`), or just manually enter for each patient.
`MakeTemplate.m` has some junk numbers that give a junk result, but should run!

"Models” contains the ODEs defining the model.
“Functions” has a few core functions which fit parameters etc. Almost every function follows a read-modify-write process, where it takes a patient struct and adds some new fields to it.
"Libraries" contains helper functions which can be mostly ignored.