function [z,u,ddu]=visualise_u()

[Re,~]=fix_all_parameters();

z=0:0.01:1;

u=Re*z.*(1-z);
ddu=-2*Re+0*z;


end


