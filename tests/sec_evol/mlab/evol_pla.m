clear all;close all;
load evolpla;
load evolpla2;

figure(1);
subplot(2,2,1)
plot(evolpla(:,1),evolpla(:,2));
xlabel('time (YR)');ylabel('h');
subplot(2,2,2)
plot(evolpla(:,1),evolpla(:,3));
xlabel('time (YR)');ylabel('k');
subplot(2,2,3)
plot(evolpla(:,1),evolpla(:,4));
xlabel('time (YR)');ylabel('p');
subplot(2,2,4)
plot(evolpla(:,1),evolpla(:,5));
xlabel('time (YR)');ylabel('q');

figure(2);
subplot(2,2,1)
plot(evolpla2(:,1),evolpla2(:,2));
hold on;
plot(evolpla2(:,1),evolpla2(:,6),'r');
xlabel('time (MJD)');ylabel('h');
subplot(2,2,2)
plot(evolpla2(:,1),evolpla2(:,3));
hold on;
plot(evolpla2(:,1),evolpla2(:,7),'r');
xlabel('time (MJD)');ylabel('k');

subplot(2,2,3)
plot(evolpla2(:,1),evolpla2(:,4));
hold on;
plot(evolpla2(:,1),evolpla2(:,8),'r');
xlabel('time (MJD)');ylabel('p');

subplot(2,2,4)
plot(evolpla2(:,1),evolpla2(:,5));
hold on;
plot(evolpla2(:,1),evolpla2(:,9),'r');
xlabel('time (MJD)');ylabel('q');
