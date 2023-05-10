function [Mw]=calculateMoment(models)


include_dipslip = 1;

muTopRange=[2.06e10 4.13e10];
muBaseRange=[3.10e10 4.13e10];
muTopRange=[3e10 3e10];
muBaseRange=[3e10 3e10];
cycle_length = 330;
rupture_length = 380e3;
faultsliprate=27e-3;
faultdiprate=[8e-3 12e-3];
faultdiprate=[10e-3 10e-3];

muTop = rand(1,length(models))*(muTopRange(2) - muTopRange(1)) + muTopRange(1);
muBase = rand(1,length(models))*(muBaseRange(2) - muBaseRange(1)) + muBaseRange(1);
dipslip = rand(1,length(models))*(faultdiprate(2) - faultdiprate(1)) + faultdiprate(1);

basedipslip = dipslip - abs(models(9,:));

if include_dipslip == 1
    topslip = sqrt((faultsliprate.^2)+(dipslip.^2)) * cycle_length;
    baseslip = sqrt(((faultsliprate - models(8,:)).^2)+(basedipslip.^2)) * cycle_length;
else
    topslip = faultsliprate * cycle_length;
    baseslip = (faultsliprate - models(8,:)) * cycle_length;
end

Mo = muTop .* (models(3,:)./sind(-models(4,:))) .* rupture_length .* topslip + muBase .* models(2,:) .* rupture_length .* baseslip;

Mw = (2/3) * (log10(Mo) - 9.1);