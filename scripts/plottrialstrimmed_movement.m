datafile = './mabel_reaching_5-4-10.mat';
binsize = 1/100;
dt = 1/10;
unitidx = 13;
processed = preprocess_movementinit(datafile, binsize, dt, unitidx);

intrial = processed.intrial;
cursorintrial = processed.cursor(intrial==1,:);
gripintrial = processed.grip(intrial==1);

figure
%Compute auto- and cross-correlation in torque and example firing rate
samplerate = 100;
binsize = 1/samplerate;
maxlag = 20;
autotorqueX = xcov(cursorintrial(:,1),samplerate*maxlag);%, 'coeff');
autotorqueY = xcov(cursorintrial(:,2),samplerate*maxlag);%, 'coeff');
autotorqueZ = xcov(cursorintrial(:,3),samplerate*maxlag);%, 'coeff');
autotorqueGrip = xcov(gripintrial(:),samplerate*maxlag);%, 'coeff');
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
saveplot(gcf, ['stim_autocorrelations_trimmed_movement.eps'], 'eps', [6 6]);