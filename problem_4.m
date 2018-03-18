parallel=[10.433 8.448 7.794 7.691 7.476 7.472 7.487];
serial=23.465;
processors=2:8
speedup=serial./parallel
plot(processors,speedup)
xlabel('processors')
ylabel('speedup')
%%
guided=mean([.174 .179 .167 .177 .178])
static=mean([.174 .167 .178 .177 .174])
dynamic=mean([.180 .176 .157 .172 .182])
