pv = 0;
pt = 0;
tt = 0;
trt = 0';

fun = @(x) 1*cos(1*x - pi()/4) + 2*cos(0.1*x - pi()/16);
[v,~,t,~,rt] = test_asq(fun,-2*pi(),2*pi(),10^-6,30);
pv = pv + v;
pt = pt + t;
tt = tt + 1;
trt = trt + rt;

fun = @(x) x.^2;
[v,~,t,~,rt] = test_asq(fun,-10,10,10^-6,30);
pv = pv + v;
pt = pt + t;
tt = tt + 1;
trt = trt + rt;

fun = @(x) x.^2;
[v,~,t,~,rt] = test_asq(fun,-10,10,10^-6,30);
pv = pv + v;
pt = pt + t;
tt = tt + 1;
trt = trt + rt;

fun = @(x) x.^3;
[v,~,t,~,rt] = test_asq(fun,-10,10,10^-6,30);
pv = pv + v;
pt = pt + t;
tt = tt + 1;
trt = trt + rt;

fun = @(x) exp(x);
[v,~,t,~,rt] = test_asq(fun,-1000,0,10^-6,30);
pv = pv + v;
pt = pt + t;
tt = tt + 1;
trt = trt + rt;

fun = @(x) log(x);
[v,~,t,~] = test_asq(fun,.1,10,10^-6,30);
pv = pv + v;
pt = pt + t;
tt = tt + 1;
trt = trt + rt;

disp('ASQ Test Suite:') 
disp(['Tests Passed: ',sprintf('%d',pv),'/',sprintf('%d',tt)])
disp(['Performance Tests Passed: ',sprintf('%d',pt),'/',sprintf('%d',tt)])
disp(['Average Performance Increase: ',sprintf('%.4f',trt/tt)])
