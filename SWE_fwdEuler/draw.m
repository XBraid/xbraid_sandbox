numfiles = 100
nspace = 41
 X = zeros(nspace,numfiles);
for k = 0:numfiles
  myfilename = sprintf('%d', k);
  X(:,k+1) = importdata(myfilename);
end
for i = 1:numfiles
#plot(X(:,i))
#sum(X(:,i))
#axis([0 nspace 0 4]);
#drawnow
#pause(0.1);
end

shallow_water_1d;
Y = ans;
for i = 1:numfiles
#plot(ans(:,i))
#axis([0 nspace 0 4]);
#drawnow
#pause(0.1);
end