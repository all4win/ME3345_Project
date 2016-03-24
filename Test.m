A = zeros(3072, 3072);
B = zeros(3072, 1);

for i = 1 : 3072
    for j = 1 : 3072
        A(i,j) = round(rand()*10);
    end
    B(i, 1) = round(rand()*10);
end

Ai = inv(A);
C = Ai * B;