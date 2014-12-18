function [tt xx yy] = orbit(r_per,r_ap,dt,Ttot)
% ORBIT  Calculates Kepler orbit as a function of time.
%   Orbit starts at perihelion at t = 0.
%   
%   Example
%   [t,x,y] = orbit(r_per,r_ap,dt,Ttot,method);
%       Output:
%           t - time, s
%           x - x in GSE coordinates, m
%           y - y in GSE coordinates, m
%       Input:
%           r_per - perigee, m
%           r_ap - apogee, m
%           dt - timestep, s
%           Ttot - total time to calculate the orbit for. If Ttot is larger
%                  than the period of the orbit, the orbit is only
%                  duplicated. It is assumed the orbit is fixed wrt the 
%                  star system, so a precession frequency (?) is added
%                  that displaces the orbit slightly.

% Physical parameters
G = 6.67384e-11; % m^3 kg^-1 s^-2 (N m^-2 kg^-2)
ME = 5.9722e24; % kg
%RE = 6371*1e3; % m

% Orbital parameters
a = (r_ap + r_per)/2; % m, semi-major axis
mu = G*ME; % m^3 s^-2
e = (r_ap - r_per)/(r_ap + r_per); % eccentricity
n = sqrt(mu/a^3);
T = 2*pi/n; % orbital period

if Ttot<T; t = 0:dt:Ttot;   
else t = 0:dt:T;
end
nt = numel(t);
tau = 0; % time of perihelion passage
M = n*(t-tau); % mean anomaly;

funEM = @(E,M) E-e*sin(E)-M;
E = fsolve(@(E) funEM(E,M),M);
r = a*(1-e*cos(E)); % radial distance from Earth
f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2)); % true anomaly (angle from x,y=0)

if Ttot<T; 
    rr = r;
    ff = f;
    tt = t;
else % copy orbit
    nT = ceil(Ttot/T);
    rr = repmat(r,1,nT);
    ff = repmat(f,1,nT);
    tt = linspace(0,nT*T,nT*nt);
end

% add precession of orbit as Earth goes around Sun
% 2*pi rad degrees in 60*60*24*365 seconds (one year)
ffplus = tt*2*pi/(60*60*24*365);
xx = rr.*cos(ff+ffplus);
yy = rr.*sin(ff+ffplus);

% remove again that extra odd part of an orbit
xx(tt>Ttot)=[];
yy(tt>Ttot)=[];
tt(tt>Ttot)=[];