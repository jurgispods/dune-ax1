function p = generate_permutation(N)

conIndex = 1;
potIndex = N * 3/4 + 1;
count = 0;
p = 0*(1:N);
size(p)

for (i=1:N)
    if (count < 3)
      p(i) = conIndex;
      conIndex = conIndex+1;
      count = count+1;
    else 
      p(i) = potIndex;
      potIndex = potIndex+1;
      count = 0;
    end


end