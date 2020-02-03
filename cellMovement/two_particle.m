function [ux,uy] = two_particle(E,sigma)

%%%  Function to calculate the displacement field in an elastic medium caused
%%%  by two flexible spherical particles' movement. 
%%%  The steps here include:
%%% 1)make two distance maps for each particle
%%% 2)make masks for each of them
%%% 3)physical interaction; force estimate
%%% 4)Green's function and calculation of ux,uy for all points
%%% 5)all the graphs, etc


