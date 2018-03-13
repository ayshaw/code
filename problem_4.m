parallel=[10.433 8.448 7.794 7.691 7.476 7.472 7.487];
serial=23.465;
processors=2:8
speedup=serial./parallel
plot(processors,speedup)
xlabel('processors')
ylabel('speedup')
