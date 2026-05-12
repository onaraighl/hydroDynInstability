function [u_vec,ddu_vec,z]=get_u()

[Re,N_1]=fix_all_parameters();

jj=1:(N_1-3);
eta=cos((pi/(N_1-2))*jj);
z=(1/2)*(eta+1);

u_vec=Re*z.*(1-z);
ddu_vec=-2*Re+0*z;

end

%**************************************************************************
