% Data zeefanalyses omrekenen naar cumulatieve zeefkrommen en percentielen
% Maarten Kleinhans augustus 2007

close all
clear
clc

% addpath('../Data')

% invoerfiles:
%M+1 zeefdiameters in mm als kolom van hoog naar laag, alle de BOVENzeef (met onderste of onderbak erbij)
%M gewichten of percentages (niet cumulatief) in zeefvolgorde in kolommen voor N monsters

%cd 'E:\onderwijs\MScVeldwerk\ZeelandGeulen';
ZeefData = load('lapeyne.txt');
ZeefDiam = load('diamlapeyne.txt');
% ZeefData = Sieveanalysis03122020S1.MassRetained(1:end-1);
% ZeefDiam = Sieveanalysis03122020S1.SieveOpeningmu*1000;
%ZeefData = load('Zuigmonsters ABS vergelijkingsproeven in GWK.txt');
ZeefDiam = ZeefDiam./1000;
[M,N] = size(ZeefData); %zie hierboven

% fractiediameters berekenen
maxtel=length(ZeefDiam);
FracDiam = NaN(length(ZeefDiam)-1,1);
for tel=1:(maxtel-1)
    FracDiam(tel) = exp((log(ZeefDiam(tel))+log(ZeefDiam(tel+1)))/2);
end %of for Fracdiam berekenen

%fracties (niet-cumulatieve kans) (%/100) in kolommen voor diverse monsters, en totaal 1 (100%) maken
i = find(ZeefData<0);
ZeefData(i) = 0;
Factor = 1./sum(ZeefData);
for tel = 1:maxtel-1
    ZeefData(tel,:) = ZeefData(tel,:).*Factor;
end %van correctie totaal P naar 100
%ZeefDat = ZeefData(2:31,:); %voor opslaan bij bijbehorende fractiediameters

%ZeefPerc = onderschrijdingskans uitrekenen (cumulatieve zeefkrommen)
ZeefPerc = zeros(size(ZeefData));
for tel = 1:maxtel-1
    for teln = 1:N
        ZeefPerc(tel,teln) = sum(ZeefData(tel:maxtel-1,teln));
    end %of for go through N samples
end %of for cumu berekening
%ZeefPerc(:,1)

%Percen (percentielen) uitrekenen
%Perc = [0.1 0.16 0.35 0.5 0.65 0.84 0.9]'; %te berekenen percentielen
Perc = [0.1 0.16 0.5 0.84 0.9 0.95]'; %te berekenen percentielen
nzeef = size(ZeefPerc);
Percen = NaN(length(Perc),nzeef(2)); %matrix voor resultaten
for dat = 1:nzeef(2) %aantal zeefkrommen
    for per = 1:length(Perc) %aantal percentielen
        i=find(ZeefPerc(:,dat)>Perc(per));
        if length(i)<1
        end %of if check D90 niet beneden tweede zeef
        top = length(i); bot = top+1;
        diam = exp(log(ZeefDiam(top))-(ZeefPerc(top,dat)-Perc(per))./...
        (ZeefPerc(top,dat)-ZeefPerc(bot,dat)).*(log(ZeefDiam(top))-log(ZeefDiam(bot))));
        Percen(per,dat) = diam;
    end %of for aantal percentielen
end %of for aantal zeefkrommen

%Moments: gemiddelde diameter en echte standaarddeviatie
GeoMeanStd= NaN(2,nzeef(2)); %matrix voor resultaten
PsiMeanStd = GeoMeanStd;
psi = log(FracDiam)./log(2); %arithmetic grain size waarbij sedimentologische phi=-psi zodat psi oploopt met D
for dat = 1:nzeef(2) %aantal zeefkrommen
    PsiMeanStd(1,dat) = nansum(psi.*ZeefData(:,dat));
    GeoMeanStd(1,dat) = 2^PsiMeanStd(1,dat);
    PsiMeanStd(2,dat) = sqrt( nansum( ((psi-PsiMeanStd(1,dat)).^2) .*ZeefData(:,dat) ) );
    GeoMeanStd(2,dat) = 2^PsiMeanStd(2,dat);
end %of for aantal zeefkrommen



%Bimodality: percentage grind en B=0.5min(Xg,Xs)/(1-Xg-Xs) met Xg=fraction>4mm Xs<1mm
% Bimod = NaN(2,nzeef(2)); %matrix voor resultaten
% igravel2mm = find(ZeefDiam(1:end-1)>=2);
% igravel4mm = find(ZeefDiam(1:end-1)>=4);
% isand1mm = find(ZeefDiam(1:end-1)<=1);
% for dat = 1:nzeef(2) %aantal zeefkrommen
    % Bimod(1,dat) = ZeefPerc(igravel2mm(end),dat);
    % Bimod(2,dat) = 0.5*(nanmin(ZeefPerc(igravel4mm(end),dat),ZeefPerc(isand1mm(1),dat))) ./...
    % (1 -ZeefPerc(igravel4mm(end),dat) -ZeefPerc(isand1mm(1),dat) ); %als in Frings-Kleinhans
% end %of for aantal zeefkrommen

%voorbeeld voor OPSLAAN
%save fracdiamExtra.txt FracDiam -ASCII;
%save fracExtra.txt ZeefData -ASCII;
%save zeefExtra.txt ZeefPerc -ASCII;
%save percenExtra.txt Percen -ASCII;
% save modder PsiMeanStd GeoMeanStd FracDiam ZeefDiam ZeefData ZeefPerc Percen Perc

%FIGURES
figure %cumulatieve zeefkromme met percentielen en gemiddelde
semilogx(ZeefDiam(1:end-1),ZeefPerc(:,1))
hold on
plot(Percen(:,1),Perc,'.')
plot(GeoMeanStd(1,1),0.5,'o',... geometric mean grain size
2^(PsiMeanStd(1,1)-PsiMeanStd(2,1)),0.5,'>',... minus one standard deviation
2^(PsiMeanStd(1,1)+PsiMeanStd(2,1)),0.5,'<')% plus one standard deviation
xlabel('particle diameter (mm)')
ylabel('Probability')

figure %niet-cumulatieve kansverdeling
semilogx(FracDiam,ZeefData(:,1))
hold on
plot(GeoMeanStd(1,1),0.01,'o',... geometric mean grain size
2^(PsiMeanStd(1,1)-PsiMeanStd(2,1)),0.01,'>',... minus one standard deviation
2^(PsiMeanStd(1,1)+PsiMeanStd(2,1)),0.01,'<')% plus one standard deviation
xlabel('particle diameter (mm)')
ylabel('fraction')
