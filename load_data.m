function [A,N,eld_y,scal_press,pitot,t_atm,M,aoa,p_atm,rho_atm,aux_inf,...
    eld_x,v_inf,v  ] = load_data( filename )

data=load(filename);
p_atm=data(:,1);			%Atmospheric pressure (Pa)
t_atm=data(:,2);			%Atmospheric Temperature (K)
rho_atm=data(:,3);			%Atmospheric Density (kg/m^3)
v_inf=data(:,4);			%Airspeed (m/s)
pitot=data(:,5);			%pitot dynamic pressure (Pa)
aux_inf=data(:,6);			%Aux dynamic pressure (Pa)
scal_press=data(:,7:21);	%scalnvale pressure
aoa=data(:,23);				%angle of attack (deg)
N=data(:,24);				%sting normal force (N)
A=data(:,25);				%sting axial force (N)
M=data(:,26);				%sting moment (Nm)
eld_x=data(:,27);			%eld probe x axis (mm)
eld_y=data(:,28);			%eld probe y axis (mm)
%v=sqrt((2.*abs((pitot-p_atm)))./rho_atm);%wake velocity (m/s)
v = sqrt((2.*aux_inf) ./ rho_atm );

end

