addpath('C:\Users\ebras\Documents\MATLAB\Add-Ons\Toolboxes\ustb');
datapath = 'D:\Master\Beamforming\TilEmil';
filename_out = 'TilEmil_'; % Name of uff file. NB different transmits are stored separately, numbered starting with 1
name = 'ch_data'; % Name of "property" 
%% Example, reading uff channel data and beamform with multiple tx/rx combinations
N = 2;
for n = 1:N
    channel_data(n) = uff.channel_data();
    channel_data(n).read(fullfile(datapath,[filename_out num2str(n) '.uff']),['/' name])
end
 
xStart = -( (channel_data(1).probe.N-1)/2 )*channel_data(1).probe.pitch ;
xEnd = xStart + (channel_data(1).probe.N-1)*channel_data(1).probe.pitch;

zStart = 0.005;
zEnd = 0.025;

lambda = channel_data(1).sound_speed/channel_data(1).pulse.center_frequency;         

dx = lambda/2;
dz = lambda/2;

bfPars.x_axis=linspace( xStart, xEnd, ceil( ((xEnd-xStart)/dx)/2)*2 ).';               % x vector [m]
bfPars.z_axis=linspace( zStart, zEnd , ceil( ((zEnd-zStart)/dz)/2)*2).';              % z vector [m] 

sca=uff.linear_scan('x_axis',bfPars.x_axis, 'z_axis', bfPars.z_axis);
   
%% Beamforming
for n = 1:N
    txangles_t(n) = channel_data(n).sequence.source.azimuth;
end

txangles = [txangles_t(1) txangles_t(1) txangles_t(1) txangles_t(1) txangles_t(1) txangles_t(2) txangles_t(2) txangles_t(2) txangles_t(2)];
rxangles = [-15 -3.5 6 15 0 -6 3.5 15 0]/15*txangles_t(2); 
txRx{1} =[-15 -3.5 6 15 0]/15*txangles_t(2); 
txRx{2} = [-6 3.5 15 0]/15*txangles_t(2);% 

bfctr = 1;
LRdata = zeros(sca.N_z_axis,sca.N_x_axis,length(txRx{1})+length(txRx{2}),size(channel_data(1).data,4));

for n = 1:N
    bfVec = txRx{n};     
    bmf_multi=pipeline();

    bmf_multi.channel_data=channel_data(n);
    bmf_multi.scan=sca;

    for bn = 1:length(bfVec)
        bmf_multi.transmit_apodization.window=uff.window.none; % Use each low resolution image separately --> no "synthetic" apod on transmit
        bmf_multi.receive_apodization.tilt = [bfVec(bn) 0];

        if bfVec(bn) == 0           
            bmf_multi.receive_apodization.window=uff.window.('hamming');
            bmf_multi.receive_apodization.f_number=1.1;
        else
            bmf_multi.receive_apodization.window=uff.window.hamming;
            bmf_multi.receive_apodization.f_number=1.4;
        end

        bmf_multi.receive_apodization.focus = sca;

        das=midprocess.das();
        das.code = code.mex;
        b_data_multi=bmf_multi.go({das});

        LRdata(:,:,bfctr,:) = reshape(b_data_multi.data,[sca.N_z_axis,sca.N_x_axis,1,size(channel_data(1).data,4)]);

        bfctr = bfctr+1;
    end
end 

% Some more useful parameters
bfPars.PRF = channel_data(1).PRF;
bfPars.c = channel_data(1).sound_speed;
bfPars.f_demod = channel_data(1).modulation_frequency;

% Speckle tracking data
preCData = LRdata(:,:,(rxangles==0),:);
compdata = sum( preCData, 3);

% Doppler data
remInds = [5 9]; 
txangles(remInds) = [];
rxangles(remInds) = [];
LRdata(:,:,remInds,:) = [];


%% PW Doppler
viewAngle = 1;

useVerasonicsdata =0;
useGEdata = 1;

addpath('C:\Users\ingvilek\vectorflow\Code\DopplerFunctions\') 

figure(); 
subplot(2,1,1);
imagesc(10*log10(abs(compdata(:,:,1)))); title('Select PW point')
pwcoord = round(ginput(1));

PWpars = bfPars;
PWpars.typeFilter = 'FIR138_3';

PWsig = squeeze(LRdata(pwcoord(2),pwcoord(1),viewAngle,:));
PWsig = hp(PWsig,PWpars);

[PW,vaxis,taxis] = PWspec(PWsig,PWpars);

subplot(2,1,2);imagesc(taxis,vaxis,10*log10(abs(PW)));
title('Select frame')
framecoords = ginput(1);

ix = floor(framecoords(1)*PWpars.PRF);
    