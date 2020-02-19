function [sys] = PtDataRead(sys)
%gets all the reqd patient data

%import raw .xls data from the patient 
[num, txt] = get_patient_data(sys.Data.PtNo);

%read times from the trial
[sys.trial_start_t, sys.trial_end_t] = get_trial_data(txt);
sys.trial_max_t = minutes(sys.trial_end_t - sys.trial_start_t);

%read cpeptide data and set the simulation time
[sys.Data.Cpep_time, sys.Data.Cpep, sys.sim_start_t, sys.sim_end_t] = get_cpeptide_data(num, txt);
sys.sim_max_t = minutes(sys.sim_end_t - sys.sim_start_t);

%read insulin bolus time and size. also convert from U to mU
[bolus_time, bolus_size] = get_bolus_data(num, txt);
sys.SC.T = minutes(bolus_time - sys.sim_start_t);
sys.SC.Ibolus = bolus_size * 1000; % convert to mU

%reading in meal times etcccc
[sys] = get_meal_data(num, txt, sys);

%Read cgm data files
[sys.Data.bg1_time, sys.Data.bg1, sys.Data.bg2_time, sys.Data.bg2] = get_cgm_data(sys.Data.PtNo);

%Read blood glucose (fingerprick) data
[sys.Data.bg3_time, sys.Data.bg3] = get_bgl_data(num, txt);

%read plassma insulin data
[sys.Data.PlasmaI_time, sys.Data.PlasmaI] = get_insulin_data(num, txt);

%Read Fasting Blood Glucose
 % fasting blood glucose level (mmol/litre)
sys.GC.fasting_bg1 = num(11, 2);
sys.GC.fasting_bg2 = num(12, 2);
sys.Data.pt_mass = num(4,2);

end