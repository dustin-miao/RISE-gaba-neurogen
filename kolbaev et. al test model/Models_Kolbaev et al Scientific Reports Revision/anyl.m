clear all 

offcl = importdata('c0.dat');
loc0 = importdata('loc1.dat');
offcl1 = importdata('c1.dat');
offcl2 = importdata('c3.dat');
offcl3 = importdata('c4.dat');

%t = 0:0.025:1000;
dt = 0.25
t100 = (100/dt)+1;
t500 = (500/dt)+1;
t800 = (800/dt)+1;
t1000 = (1000/dt)+1;

plot(loc0,offcl(:,t100),':','Color','red')
hold on 
plot(loc0,offcl(:,t500),':','Color','black')
hold on 
plot(loc0,offcl(:,t800),':','Color','blue')
hold on 
plot(loc0,offcl(:,t1000),':','Color','green')
box off 
hold on 

plot(loc0,offcl1(:,t100),'red')
hold on 
plot(loc0,offcl1(:,t500),'black')
hold on 
plot(loc0,offcl1(:,t800),'blue')
hold on 
plot(loc0,offcl1(:,t1000),'green')
box off 

plot(loc0,offcl2(:,t100),'--','Color','red')
hold on 
plot(loc0,offcl2(:,t500),'--','Color','black')
hold on 
plot(loc0,offcl2(:,t800),'--','Color','blue')
hold on 
plot(loc0,offcl2(:,t1000),'--','Color','green')
box off 
hold on 

plot(loc0,offcl3(:,t100),'-o','Color','red')
hold on 
plot(loc0,offcl3(:,t500),'-o','Color','black')
hold on 
plot(loc0,offcl3(:,t800),'-o','Color','blue')
hold on 
plot(loc0,offcl3(:,t1000),'-o','Color','green')
box off 
