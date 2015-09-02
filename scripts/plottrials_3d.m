load('mabel_reaching_5-4-10')

x = Cursor_X;
y = Cursor_Y;
z = Cursor_Z;

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

%Whole recording lasts 6000s...
tstart = 555.64; tend = 570.64;
tidx = (tstart*100):(tend*100);
tt = tidx/100;

eventbins = zeros(size(Cursor_X));

%Plot Events
%Only plots events within range
idxstart = find(Events_Data(1,:)>=(tstart*1000), 1);
idxend = find(Events_Data(1,:)>=(tend*1000), 1);
evtbins = ceil(Events_Data(1,idxstart:idxend)/10);
evts = Events_Data(2,idxstart:idxend);
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
	evttime = evtbins(idx);
	if evt == ORIGINAPPEARS;
		inorigin = 1;
		originstart = evttime;
	elseif evt == ORIGINREACHED;
		if ~inorigin
			originstart = tidx(1);
		end
		eventbins(originstart:evttime) = 1;
		inorigin = 0;
	elseif evt == TARGETAPPEARS;
		intarget = 1;
		targetstart = evttime;
	elseif evt == TARGETREACHED;
		if ~intarget
			targetstart = tidx(1);
		end
		eventbins(targetstart:evttime) = 2;
		intarget = 0;
	elseif evt == GOSIGNAL;
		ingo = 1;
		gostart = evttime;
	elseif evt == RELEASEGRIP;
		if ~ingo
			gostart = tidx(1);
		end
		eventbins(gostart:evttime) = 3;
		ingo = 0;
	end
end 

if inorigin
	eventbins(originstart:tidx(end)) = 1;
end
if intarget
	eventbins(targetstart:tidx(end)) = 2;
end
if ingo
	eventbins(gostart:tidx(end)) = 3;
end


%Plot cursor motion in 3D
figure
hold on 
scatter3(x(tidx), y(tidx), z(tidx), 100*Grip_force(tidx), eventbins(tidx));
plot3(x(tidx), y(tidx), z(tidx), 'k');


saveplot(gcf, './trialsample_3d.eps', 'eps', [6 6])
