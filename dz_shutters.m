% Shutter EPI excitations.
% Assumes 90 degree tip; 
% lower flips should switch to 'st' design of rfShut, 
% and scale flip accordingly.

% total # of shutters is R*Nshots
% so the width of each shutter should be FOV/(R*Nshots),
% and the phase encode XFOV should be (Nshots)*(shutter width)

% rfSl: The slice-selective pulse (radians; scaled to 1)
% rfShut: The shutter envelope (radians; scaled to pi/2)
% gpos: The slice-select trapezoid (g/cm)
% gyBlip: The blip between the subpulses (g/cm)
% rfEP: The total constructed pulse (radians)
% gzRew: the slice-select rewinder
% gyRew: the phase encode/slab-select rewinder
% ttipdown: the time into the pulse when the TE period should start
% rfFM: normalized FM waveform for slice shifting
% rfPhs: matrix of RF phase waveforms that shift the shutter for each shot
% phsMtx: matrix of RF phases (one per subpulse) that shift the shutter for
%         each shot

% dependencies:
% - JP's RF toolbox (for dzrf)
% - domintrap.m - design trapezoid for a target area under the flat
%   portion (for non-ramp-sampled pulses).
% - dotrap.m - design trapezoid for total target area (for rewinders)
% - Fessler IRT for im (could replace with imagesc)
% - rf2ppe.m - write rf pulse to Philips PPE
% - gr2ppe.m - write gradient pulses to Philips PPE

dt = 6.4e-6; % s, dwell times of RF + gradient pulses
Nshots = 3; % number of shots/EPI segments
imFOV = 20.2; % cm, imaging FOV in shuttered dim
R = 4; % imaging acc factor
doSim = true; % do final Bloch simulation
flip = 30; % flip angle

tbw = [4 6]; % time bandwidth in slice, shutter dims
dthick = [0.2 imFOV/(R*Nshots)]; % slice thickness, shutter width

kw = tbw ./ dthick; % width of k-space coverage in each dimension (1/cm)

gz_area = kw(1) / 4257; % z(slice)-gradient area (g-s/cm)
gzmax = 4; % g/cm
gymax = 4; % g/cm
gslew = 20000; % g/cm/s
[gpos,ramppts] = domintrap(gz_area,gzmax,gslew,dt); % plateau sums to desired area
% remove last point since it is zero and will give two consecutive zeros in
% total waveform
gpos = gpos(1:end-1);
Ntz = length(gpos);

% design slice-selective subpulse
rfSl = real(dzrf(Ntz-2*ramppts+1,tbw(1),'st','ls',0.01,0.01)); % arb units
% normalize to one radian flip
rfSl = rfSl./sum(rfSl);

% design the shutter envelope
if flip == 90
    rfShut = real(dzrf(round(kw(2)*Nshots*dthick(2)),tbw(2),'ex','ls',0.01,0.01)); % radians
else
    rfShut = real(dzrf(round(kw(2)*Nshots*dthick(2)),tbw(2),'st','ls',0.01,0.01)); % arb units
    % normalize to target flip
    rfShut = rfShut./sum(rfShut)*flip*pi/180;
    %rfShut(1:2:end) = 0;
end

% construct the pulse with gaps for ramps
rfEP = kron(rfShut,[zeros(1,ramppts) rfSl zeros(1,ramppts-1)]);
ttipdown = length(rfEP)/2*dt*1000; % time into the pulse at which TE should start (ms) - calculate before we add rewinder zeros

% build total gz gradient waveform
gzEP = kron(ones(1,floor(length(rfShut)/2)),[gpos -gpos]);
if rem(length(rfShut),2)
    gzEP = [gzEP gpos];
end

% get the gy blips
gyBlip = dotrap(1/(Nshots*dthick(2))/4257,gymax,gslew,dt);
if rem(length(gyBlip),2)
    gyBlip = [gyBlip 0]; % add a zero to make it even length
end
% append zeros so that they straddle consecutive gz traps
gyBlipPad = [zeros(1,Ntz-length(gyBlip)) gyBlip];
gyEP = [zeros(1,length(gyBlip)/2) kron(ones(1,length(rfShut)-1),gyBlipPad)];
gyEP = [gyEP zeros(1,length(gzEP)-length(gyEP))];

% calculate and add rewinders
gzRew = dotrap(sum(gpos)*dt/2,gzmax,gslew,dt);
gzEP = [gzEP ((-1)^rem(length(rfShut),2))*gzRew];
gyRew = -dotrap(sum(gyBlip)*dt*(length(rfShut)-1)/2,gymax,gslew,dt);
gyEP = [gyEP gyRew];

% zero pad waveforms to same length
gzEP = [gzEP zeros(1,max(length(gzEP),length(gyEP))-length(gzEP))];
gyEP = [gyEP zeros(1,max(length(gzEP),length(gyEP))-length(gyEP))];
gEP = [gzEP(:) gyEP(:)]; % stick them together into matrix
rfEP = [rfEP(:); zeros(size(gEP,1)-length(rfEP),1)];

% calculate FM waveform for slice-shifting
rfFM = repmat([ones(Ntz,1);-ones(Ntz,1)],[floor(length(rfShut)/2) 1]);
if rem(length(rfShut),2)
    rfFM = [rfFM;ones(Ntz,1)];
end
rfFM = [rfFM; zeros(length(rfEP)-length(rfFM),1)];

% calculate the phases to shift the slab to the other locations
phsMtx = angle(exp(1i*2*pi*(0:Nshots-1)'/Nshots*(0:length(rfShut)-1)));
rfPhs = kron(phsMtx,ones(1,Ntz));
rfPhs = [rfPhs zeros(Nshots,length(rfEP)-size(rfPhs,2))];

% confirm that the shutters go where we expect
figure;hold on
y = (-64:63)/128*Nshots*dthick(2);
for ii = 1:Nshots
    plot(y,abs(fftshift(fft(rfShut.*exp(1i*phsMtx(ii,:)),128))));
    ylabel 'approx flip angle'
    xlabel 'cm'
end
title 'All shutters'

% plot the pulses
figure
subplot(411)
t = (0:length(rfEP)-1)*dt*1000; % time in ms
plot(t,rfEP);
xlabel 'ms',ylabel 'radians'
title 'RF pulse'
subplot(412)
plot(t,gEP);
xlabel 'ms',ylabel 'g/cm'
title 'Gradient waveforms'
legend('gz','gy');
subplot(413)
plot(t,rfFM)
xlabel 'ms',ylabel 'arb units'
title 'Normalized FM waveform for slice shifting'
subplot(414)
plot(t,rfPhs)
xlabel 'ms',ylabel 'Radians'
title 'Phase waveforms to shift shutter to each location (one per shot)'
c = axis; axis([c(1) c(2) -4 4]);

if doSim
    disp 'Bloch-simulating final pulses'
    mxy = blochsim_spinor(rfEP/(2*pi*4257*dt),gEP,...
        [10*dthick(1) dthick(2)*Nshots],[128 128],zeros(128),dt).';
    figure;im((-64:63)/128*dthick(2)*Nshots,(-64:63)/128*10*dthick(1),mxy);
    xlabel 'y (shutter), cm'
    ylabel 'z (slice), cm'
    mxyInPlane = zeros(101,101,Nshots); % to match XY's sims
    for ii = 1:Nshots
        mxyInPlane(:,:,ii) = blochsim_spinor(rfEP.*exp(1i*rfPhs(ii,:)')/(2*pi*4257*dt),[0*gEP(:,1) gEP(:,2)],...
            [20.2 20.2],[101 101],zeros(101),dt).';
    end
    save(sprintf('dz_shutters_tb%d_nshots%d_R%d',tbw(2),Nshots,R),'mxyInPlane');
end

% write to philips PPE
rfAng = sum(rfEP)*180/pi; % convert to degrees
rf2ppe('shutter_rf.txt',length(rfEP)*dt*1000,rfAng,ttipdown,0,rfEP,rfFM);
gr2ppe('shutter_grads.txt',size(gEP,1)*dt*1000,ttipdown,10*[gEP zeros(size(gEP,1),1)]);


