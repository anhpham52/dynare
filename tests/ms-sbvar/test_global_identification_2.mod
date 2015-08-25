// same as test_lower_cholesky.mod, but using exclusion syntax
var R Pie Y;

varobs Y Pie R;

svar_identification;
exclusion lag 0;
equation 1, Pie, Y;
exclusion lag 1;
equation 2, Y;
end;

svar_global_identification_check;

