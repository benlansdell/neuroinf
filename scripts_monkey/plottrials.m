load('mabel_reaching_5-4-10')

%Legend:
TRIALSTART = 10;
TRIALEND = 15;
ORIGINAPPEARS = 100;
ORIGINREACHED = 101;
TARGETAPPEARS = 110;
TARGETREACHED = 111;
GOSIGNAL = 103;
RELEASEGRIP = 104;
GRIPPRESSED = 368;
GRIPRELEASED = 369;

figure
clf
%Whole recording lasts 6000s...
tstart = 525.64; tend = 575.64;
tidx = (tstart*100):(tend*100);
tt = tidx/100;

%Plot grip force
hold on
h(1) = plot(tt, Grip_force(tidx), '--');
%h(1) = plot(tt, Cursor_X(tidx), '--', tt, Cursor_Y(tidx), '--', tt, Cursor_Z(tidx), '--');
colors = colormap

%Plot Events
%Only plots events within range
idxstart = find(Events_Data(1,:)>=(tstart*1000), 1)
idxend = find(Events_Data(1,:)>=(tend*1000), 1)
evttimes = Events_Data(1,idxstart:idxend)/1000;
evts = Events_Data(2,idxstart:idxend);
h(2) = plot(evttimes, ones(length(evttimes))*.01, 'x');
xlim([min(tt), max(tt)]);

%Add trial start/end lines
intrial = 0;
inorigin = 0;
intarget = 0;
ingo = 0;
ingrip = 0;
trialstart = 0;
originstart = 0;
targetstart = 0;
gostart = 0;
gripstart = 0;

for idx = 1:length(evts)
	evt = evts(idx);
	evttime = evttimes(idx);
	if evt == TRIALSTART;
		intrial = 1;
		trialstart = evttime;
	elseif evt == TRIALEND;
		if ~intrial
			trialstart = tstart;
		end
		%Draw line between here and trialstart
		h(3) = plot([trialstart, evttime], [-0.01, -0.01], 'r', 'LineWidth', 2);
		%line([trialstart, evttime], [-0.01, -0.01], 'LineWidth', 2)
		intrial = 0;
	elseif evt == ORIGINAPPEARS;
		inorigin = 1;
		originstart = evttime;
	elseif evt == ORIGINREACHED;
		if ~inorigin
			originstart = tstart;
		end
		%Draw line between here and trialstart
		h(4) = plot([originstart, evttime], [-0.03, -0.03], 'b', 'LineWidth', 2);
		%line([originstart, evttime], [-0.03, -0.03], 'LineWidth', 2)
		inorigin = 0;
	elseif evt == TARGETAPPEARS;
		intarget = 1;
		targetstart = evttime;
	elseif evt == TARGETREACHED;
		if ~intarget
			targetstart = tstart;
		end
		%Draw line between here and trialstart
		h(5) = plot([targetstart, evttime], [-0.05, -0.05], 'g', 'LineWidth', 2);
		%line([targetstart, evttime], [-0.05, -0.05], 'LineWidth', 2)
		intarget = 0;
	elseif evt == GOSIGNAL;
		ingo = 1;
		gostart = evttime;
	elseif evt == RELEASEGRIP;
		if ~ingo
			gostart = tstart;
		end
		%Draw line between here and trialstart
		h(6) = plot([gostart, evttime], [-0.07, -0.07], 'k', 'LineWidth', 2);
		%line([gostart, evttime], [-0.07, -0.07], 'LineWidth', 2)
		ingo = 0;
	elseif evt == GRIPPRESSED;
		ingrip = 1;
		gripstart = evttime;
	elseif evt == GRIPRELEASED;
		if ~ingrip
			gripstart = tstart;
		end
		%Draw line between here and trialstart
		h(7) = plot([gripstart, evttime], [-0.09, -0.09], 'y', 'LineWidth', 2);
		%line([gripstart, evttime], [-0.09, -0.09], 'LineWidth', 2)
		ingrip = 0;
	end
end 

%Close off any events unfinished:
if intrial
	h(8) = plot([trialstart, tend], [-0.01, -0.01], 'r', 'LineWidth', 2);
	%line([trialstart, tend], [-0.01, -0.01], 'LineWidth', 2)
end
if inorigin
	h(9) = plot([originstart, tend], [-0.03, -0.03], 'b', 'LineWidth', 2);
	%line([originstart, tend], [-0.03, -0.03], 'LineWidth', 2)
end
if intarget
	h(10) = plot([targetstart, tend], [-0.05, -0.05], 'g', 'LineWidth', 2);
	%line([targetstart, tend], [-0.05, -0.05], 'LineWidth', 2)
end
if ingo
	h(11) = plot([gostart, tend], [-0.07, -0.07], 'k', 'LineWidth', 2);
	%line([gostart, tend], [-0.07, -0.07], 'LineWidth', 2)
end
if ingrip
	h(12) = plot([gripstart, tend], [-0.09, -0.09], 'y', 'LineWidth', 2);
	%line([gripstart, tend], [-0.09, -0.09], 'LineWidth', 2)
end

hl = legend([h([1 3:7])], 'Grip force', 'Within trial', 'Origin seek', 'Target seek', 'Go signal', 'Grip pressed')
saveplot(gcf, './trialsample.eps', 'eps', [6 6])

%Plot autocorrelation too...
%Plot a bunch of preprocessing diagnostics
figure
%Compute auto- and cross-correlation in torque and example firing rate
samplerate = 100;
binsize = 1/samplerate;
maxlag = 90;
autotorqueX = xcov(Cursor_X(:),samplerate*maxlag);%, 'coeff');
autotorqueY = xcov(Cursor_Y(:),samplerate*maxlag);%, 'coeff');
autotorqueZ = xcov(Cursor_Z(:),samplerate*maxlag);%, 'coeff');
autotorqueGrip = xcov(Grip_force(:),samplerate*maxlag);%, 'coeff');
tt = -maxlag:binsize:maxlag;
subplot(2,2,1)
plot(tt, autotorqueX);
title('Cursor X');		
subplot(2,2,2)
plot(tt, autotorqueY)
title('Cursor Y');
subplot(2,2,3)
plot(tt, autotorqueZ);
title('Cursor Z')
subplot(2,2,4)
plot(tt, autotorqueGrip);
title('Grip force');		
saveplot(gcf, ['stim_autocorrelations.eps'], 'eps', [6 6]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Double check processed structure labels in trial times correctly%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datafile = './data/mabel_reaching_5-4-10.mat';
nU = 45;
binsize = 1/100;
samplerate = 1/binsize;
nK_sp = 40;
nK_stm = 6;
a = 10;
dt_sp = binsize;
dt_stm = 100/1000;
const = 'on';
percentage = 80;
models = {};
nRep = 20;
fn_out = './monkeyresults/run_glm_fitting_sp_100Hz_fwdrev_trim_80pt.eps';
processed = preprocess_monkey(datafile, binsize, 1);
nB = size(processed.cursor,1);

istart = 1;
iend = 5000;
ii = istart:iend;
tt = (ii)*binsize;

figure 
clf
intrial = processed.intrial(ii)==1;
plot(tt, processed.grip(ii), '--')%, tt, processed.cursor(ii, 1), tt, processed.cursor(ii, 2), tt, processed.cursor(ii, 3));
hold on
%plot(tt(intrial), processed.grip(intrial), '.', tt(intrial), processed.cursor(intrial, 1), '.',...
% tt(intrial), processed.cursor(intrial, 2), '.', tt(intrial), processed.cursor(intrial, 3), '.')

plot(tt(intrial), processed.grip(intrial), 'r.')
%Compare to images made above... seems to check out