function [f] = greenwood(n,species)
%greenwood Greenwood function for cochlea frequency mapping according to
% n: [0,1] relative cochlear length measured from the base
% f: in kHz
% Greenwood functions are in the forms:
%  F = A* (10 .^ (a*x) - k)
% where x can be expressed in physical distance (e.g. mm) or proportion
% (e.g. 0-1 or 0-100%) of basilar length, from the basal end.

if nargin < 2; species = "human";end

switch species
    case {"cat","Cat"} % ref: Greenwood (1990); Liberman (1982)
        n = 100 .* n;
        f = ( 10 .^ ( (100.0 - n) * 0.021)  - 0.8) .* 0.456;
    case {"human","Human"} % ref: ??
        n = 100 .* n;
        f = ( 10 .^ ((100.0 - n) / 100.0 * 2.0) - 0.4) * 200.0 / 1000.0;
    case {"human-greenwood1990","Human-greenwood1990"} % ref: Greenwood 1990
        A = 165.4;
        a = 2.1;
        x = 1 - n;
        k = 0.88;
        f = A .* ( 10 .^ (a * x) - k) ./ 1000;
    case {"mouse","Mouse"} % ref: modified from Müller (2004) using data from Taberner & Liberman (2005) 
        A = 9.8;
        a = 0.92;
        x = 1 - n;
        k = 0.68;
        f = A .* ( 10 .^ (a * x) - k) ;
    case {"mouse-Mueller2005","Mouse-Mueller2005"} % LINEAR!! 
        % d = 156.5 - 82.5 * log(f) (d in % length from base; f in kHz)
        d = 100 .* n;
        f = 10 .^ ( (156.6 - d) / 82.5 );
    case {"mouse-Ehret1975","Mouse-Ehret1975"} % Ehret (1975)
        A = 3.350;
        a = 0.21;
        x = 7* (1 - n); % mm
        k = 0;
        f = A .* ( 10 .^ (a * x) - k) ;
    case {"mouse-Muniak2013","Mouse-Muniak2013"} % Muniak (2013)
        % d = 78.43% * log10(f) - 49.96% % f in kHz; d in % from APEX
        % d = B * log10(f) - b
        % f = 10^(-b/B) * 10^(1/B * d)
        A = 10 ^ (49.96 / 78.43);
        a = 1/78.43;
        d = 100*(1-n);
        k = 0;
        f = A .* ( 10 .^ (a * d) - k) ;
    otherwise
        f = nan; %A .* ( 10 .^ (a * n) - k);
end

end


%Greenwood function for cochlea frequency mapping according to
% imageJ plug-in from the Eaton-Peabody Laboratories
% https://masseyeandear.org/research/otolaryngology/eaton-peabody-laboratories/histology-core
% https://meeeplfiles.partners.org/Measure_line.class
% https://masseyeandear.org/assets/MEE/pdfs/ent/ReadMe.v5.doc

% public double fCat(double n) {
%         n *= 100.0;
%         return (Math.pow(10.0, (100.0 - n) * 0.021) - 0.8) * 0.456;
%     }
%     
%     public double fGuinea_pig(double n) {
%         n *= 100.0;
%         return Math.pow(10.0, (66.4 - n) / 38.2);
%     }
%     
%     public double fChinchilla(double n) {
%         n *= 100.0;
%         return Math.exp((100.0 - n) * 0.051) * 125.0 / 1000.0;
%     }
%     
%     public double fHuman(double n) {
%         n *= 100.0;
%         return (Math.pow(10.0, (100.0 - n) / 100.0 * 2.0) - 0.4) * 200.0 / 1000.0;
%     }
%     
%     public double fMouse(final double n) {
%         return (Math.pow(10.0, (1.0 - n) * 0.92) - 0.68) * 9.8;
%     }
%     
%     public double fRat(final double n) {
%         return Math.log((n * 100.0 + 4.632) / 102.048) * -1.0 / 0.04357;
%     }
%     
%     public double fRhesusMonkey(final double n) {
%         return (Math.pow(10.0, (1.0 - n) * 2.1) - 0.85) * 360.0 / 1000.0;
%     }
%     
%     public double fGerbil(final double n) {
%         return 0.398 * (Math.pow(10.0, (1.0 - n) * 2.2) - 0.631);
%     }
%     
%     public double fMarmoset(final double n) {
%         return 255.7 * (Math.pow(10.0, (1.0 - n) * 2.1) - 0.686) / 1000.0;
%     }