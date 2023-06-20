%compare_vep_psych

clear all; close all

vp = readtable('C:\Users\Ione Fine\Documents\code\dynamic_contrast\fitdata_vep_psychophysics\vep_psychophysics_congruent_fits.csv');
v = readtable('C:\Users\Ione Fine\Documents\code\dynamic_contrast\fitdata_vep\vep_congruent_fits.csv');
group = {'AM', 'NS', 'BD'};
CData =[ 1 0 0 ; 0 0 1;  1 0 1];
vp.kAE(vp.kAE<0) = NaN;
v.kAE(v.kAE<0) = NaN;
figure(1); clf
for g = 1:3
    ind  = find(strcmp(vp.group, group{g}));
    length(ind)
    p= scatter(v.kAE(ind), vp.kAE(ind), 'o','filled'); hold on
    p.SizeData = 50;
    p.MarkerEdgeColor= CData(g,:);
    p.CData =CData(g,:);
end
set(gca, 'XLim', [0 1]);
set(gca, 'YLim', [0 1]);
plot([ 0 1], [0 1], 'k--')
av.kAE = v.kAE(ind); avp.kAE = vp.kAE(ind); 
gv =find( ~isnan(av.kAE) .* ~isnan(avp.kAE));
corr(av.kAE(gv), avp.kAE(gv))
xlabel('VEP'); ylabel('VEP Psychophysics');
axis square;

figure(2); clf
vp.uFE(vp.uFE<0) = NaN;
v.uFE(v.uFE<0) = NaN;
for g = 1:3
    ind  = find(strcmp(vp.group, group{g}));
    length(ind)
    p= scatter(v.uFE(ind), vp.uFE(ind), 'o','filled'); hold on
    p.SizeData = 50;
    p.MarkerEdgeColor= CData(g,:);
    p.CData =CData(g,:);
end
% set(gca, 'XLim', [0 1]);
% set(gca, 'YLim', [0 1]);
plot([ 0 1], [0 1], 'k--')

av.uFE = v.uFE(ind); avp.uFE = vp.uFE(ind); 
gv =find( ~isnan(av.uFE) .* ~isnan(avp.uFE));
corr(av.uFE(gv), avp.uFE(gv))
xlabel('VEP'); ylabel('VEP Psychophysics');
axis square;


