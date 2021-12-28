function [fs,success,iteration] = spaced_frequencies(separation, m)
javaaddpath java
fs = zeros(m,1);
max_iterations = 10000;
iteration = 1;
while iteration < max_iterations
  r = RandSep(separation);
  success = 1;
  for i=1:m
    fs(i) = r.next();
    if fs(i)==-1
      success = 0;
      break;
    end
  end
  if success
    break;
  end
  iteration = iteration + 1;
end
end
