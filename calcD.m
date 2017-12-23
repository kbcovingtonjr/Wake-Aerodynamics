function [dPrime] = calcD(Vinf,RHOinf,Vel,h)
%%%%%%

% Authors:      Lucas Calvert
%               Keith Covington
%               Daniel Mastick
%               Ginger Beerman

%ASEN 3111 - Experimental Lab 1
%Created: 9/28/17
%Edited: 9/28/17
%
%Inputs: Free Stream Velocity (vector), Free Stream Pressure (vector), Measured Velocities in
%the wake (vector), Corresponding heights (vector
%
%Outputs: Drag per unit span (D')
%
%This function calculates the drag per unit span for the data from one test
%%%%%%%

%First, calculate the velocity decrement at each point measured:
vDec = Vinf - Vel;

%Multiply by the measured velocity (elementwise):
v = Vel.*vDec;

%Now, numerically integrate over the height:
result = trapz(h/1000,v);   %convert to meters


%Average the free stream density data points:
RHOinf_avg = mean(RHOinf);

%Multiply by the density to get Drag per unit span:
D = RHOinf_avg*result;

%Calculating the drag coefficient per unit span
c = 0.0889; %Chord length of the airfoil [m]
dPrime = (2*D)/(RHOinf_avg*(mean(Vinf)^2)*c);
end

