function Trial = LoadData(Trial)
% patientNum: internal identifier
% patientCode: human-readable identifier

issource = @(dataset) dataset == Trial.source;

%% Select function and numbers.
nums = containers.Map;
if issource("Detemir")
    MakeDataFunction = @MakeDetemir;
    nums("all") = [1 3 4];
elseif issource("DISST")
    MakeDataFunction = @MakeDISST;
    nums("all") = [1:50];
    nums("best") = [3 5 7 8 9 13 14 16 24 25];
elseif issource("CREBRF2021")
    MakeDataFunction = @MakeCREBRF2021;
    nums("all") = [18  24  26  39  33  28  35  37  46  56  97  80 116 122 128 ...
        129 123 141 139 156 160 161 158 182 172 180 186 135 115 145 205  ...
        209 151 178 191 216 175 240 247 248 170 259  49  98  99 105   ...
        86 222   8  12  29  40  44  51  69  42  77  87  89  63 147 198 ...
        203 204 218 103 194 219  19  27 104 168 190 236 253 257 233  ...
        2   3   4   5  15  21  20  22  23   6  31  25  32  45  48   ...
        50  53  43  61  65  72  76  66  73  85  92  67  47  79 120  ...
        93 125 134 140 138 113 155 126 154 163 167 162 176 184 177  ...
        187 188 189 192  83  91 112 132 159 153 169 196 165 206 202  ...
        208 171 213 220 231 173 221 235 241 237 249 199 244 183 254 251  ...
        2   4  13  23  25  29  33  36  39  41  54  76  80  84    ...
        1   8   9  19  30  32  40  43  62  69  75   3  14  17  27  28   ...
        34  37  47  53  57  85];
    nums("best") = [1 3 4 8 9 14 17 18 19 20 23 25 27 27 28 29 30 31 32 33 34 36 37 37 40 41 42 43 45 47 48 53 54 57 69 72 75 76 76 79 83 84 85 85 86 87 91 93 99 103 105 112 120 135 140 141 145 151 153 155 156 158 160 162 167 169 171 172 173 175 177 178 180 182 184 186 188 189 190 194 198 199 203 204 205 206 213 219 220 221 222 231 233 235 236 240 241 244 247 249 251 253 257];
elseif issource("CREBRF")
    MakeDataFunction = @MakeCREBRF;
    nums("all") = [146 95 68 140 12 19 147 154 33 85 126 46 156 104 72 79 ...
        73 65 78 105 138 158 87 198 128 169 186 153 115 209 196 160 145 ...
        216 166 171 220 259 240 253 235 194 263 251];
    nums("best") = [12 128 146 160 166 169 171 196 198 216];
elseif issource("OGTTLui")
    MakeDataFunction = @MakeOGTTLui;
    nums("all") = [1 2 4 5 14 16 22 23 25 30];
    nums("best") = [1 2 4 14 22 23 25 30];
end

% Replace numbers if selecting all/best.
if istext(Trial.patients) && any(contains(keys(nums), Trial.patients))
    patientNums = nums(Trial.patients);
else
    patientNums = Trial.patients;
end

%% Load data.
patientSet = MakeDataFunction(Trial, patientNums);

%% Add remaining elements to patient structs.
for ii = 1:length(patientSet)
    patientSet{ii}.options.ODE = Trial.Config.DEFAULTODEOPTIONS;
    patientSet{ii}.figures = gobjects(0);
    
    patientSet{ii} = AddParameters(patientSet{ii});
    patientSet{ii} = AddTrialInputs(patientSet{ii});
end

%% Set up output directory.
recipeStruct = functions(Trial.recipe);
recipeName = string(recipeStruct.function);
trialPath = fullfile(Trial.source, recipeName);

% Append label if present.
if isfield(Trial, "label")
    trialPath = fullfile(trialPath, Trial.label);
else
    Trial.label = "";
end

Trial.outputPath = trialPath;

%% Save and return.
Trial.figureList = DefaultFigures();
Trial.patientSet = patientSet;
Trial.timepoint = tic;

end

