function [a15avg,a25avg,c15avg,c25avg] = compareDrag(D)
%%%%%%%
% Authors:      Lucas Calvert
%               Keith Covington
%               Daniel Mastick
%               Ginger Beerman

%ASEN 3111 - Experimental Lab 1
%Created: 10/15/17
%Edited: 10/15/17
%
%Inputs: drag per unit span for every test
%
%Outputs: Average drag per unit span for each test type 
%
%This function compares the drag coefficient per unit span for similar tests
%%%%%%%

%% Initializing to zero
a15=[];
a25=[];
c15=[];
c25=[];

%% Grouping similar tests
for(i=0:2) %Number of section
a15 = [a15 D([(13+40*i) (14+40*i) (17+40*i) (18+40*i) (21+40*i) (22+40*i)...
    (25+40*i) (26+40*i) (29+40*i) (30+40*i) (33+40*i) (34+40*i) (37+40*i)...
    (38+40*i)])];
a25 = [a25 D([(15+40*i) (16+40*i) (19+40*i) (20+40*i) (23+40*i) (24+40*i)...
    (27+40*i) (28+40*i) (31+40*i) (32+40*i) (35+40*i) (36+40*i) (39+40*i)...
    (40+40*i)])];

c15 = [c15 D([(1+40*i) (2+40*i) (5+40*i) (6+40*i) (9+40*i) (10+40*i)])];
c25 = [c15 D([(3+40*i) (4+40*i) (7+40*i) (8+40*i) (11+40*i) (12+40*i)])];
end

%% Excluding tests that are not in the data set
a15(isnan(a15)) = [];
a25(isnan(a25)) = [];
c15(isnan(c15)) = [];
c25(isnan(c25)) = [];

%% Averaging the values
a15avg = mean(a15);
a25avg = mean(a25);
c15avg = mean(c15);
c25avg = mean(c25);

end

