%% ME 3345
%% Spring 2015
%% Team project -- Topic 2 (Heat sink problem)
%% AUTHOR:  Tiancheng Gong

%% PURPOSE: 
%      Find and plot the temperature distribution of the heat sink, find
%  the position with the highest temperature using numerical method.

%% ASSUMPTIONS: 
%      - The heat transfer coefficient is uniform.
%      - The heat flux provided by the generation in the chip is uniform.
%      - Steady state.
%      - 3D convection.
%      - No radiation happens.
%      - The size of the block unit used is .5mm*.5mm*.5mm. The temperature
%      within in each block is uniform.
%      - The temperature of the heat sink is around 300K.
%      - The temperature of the infinite distance is the base, with the
%        value of 0, since the only thing we care is the temperature
%        difference between the chip and the environment. IT DOES NOT MEAN
%        THE TEMPERTURE IS 0 K.

%% IMPORTANT VALUES:
%      k = 237 W/m*K;
%      h = 100 W/m^2*K;
%      size = 0.5 mm = 5e(-4) m;
%      Total power = 4W;

%% ANALYSIS AND METHOD:
%      - Since the shape of the heat sink is symmetric, the model can be
%        simplized by dividing the heat sink into 4 identical parts, 4
%        squares and one rod in the middle of each square. (It can be
%        divided into at most 8 identical triangle parts, but it is not
%        easy to build this model using cube blocks.)
%      - The tool used here is the finite difference equations. For each
%        cube block unit, we can derive a finite difference equation based
%        on the temperatrue of its neighbors and the property of the
%        material. We will have n cubes as n unknowns as well as n
%        equations, thus the equations are solvable for any temperature
%        difference between the chip and the environment.
%      - Using trial and error, we will finally get a reasonable
%        temperature difference for the power of 4W. Using the solution
%        derived, we can thus plot the temperature distribution and find
%        the position of the heighest temperature.


%% PART 1: Create the array for the coefficients of the equations.
%  Linear Equations: At = B; B is all zeros.
size = .0005;
a = .008 / size;
b = .008 / size;
c = .006 / size;
total_cubes = a * b * c;
A = zeros(total_cubes, total_cubes);
B = zeros(total_cubes, 1);

% Set the value of coefficients based on finite difference theorem.
% Consider the border cubes carefully.
area = size ^ 2;
k = 237;
h = 100;
heat_from_chip = 1 / (a * b);
% for i = 1 : total_cubes
%     x = mod(i, 16);
%     y = floor(i / (a * b));
%     z = floor(mod(i, (a * b)) / a);
%     if (y == 0)
%         B(i, 1) = heat_from_chip;
%     end
%     if (y == (c - 1))
%     end
% end
for x = 0 : (a - 1)
    for y = 0 : (b - 1)
        for z = 0 : (c - 1)
            cur_index = z * (a * b) + y * a + x + 1;
            % left neighbor
            if (x ~= 0)
                left_index = cur_index - 1;
                coe = k / size * area;
                A(cur_index, left_index) = - coe;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
            else
                coe = h * area;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
            end
            % right neighbor
            if (x ~= (a - 1))
                right_index = cur_index + 1;
                coe = k / size * area;
                A(cur_index, right_index) = - coe;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
            end
            % front neighbor
            if (y ~= 0)
                front_index = cur_index - a;
                coe = k / size * area;
                A(cur_index, front_index) = - coe;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
            else
                coe = h * area;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
            end
            % back neighbor
            if (y ~= (b - 1))
                back_index = cur_index + a;
                coe = k / size * area;
                A(cur_index, back_index) = - coe;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
            end
            % down neighbor
            if (z ~= 0)
                down_index = cur_index - a * b;
                coe = k / size * area;
                A(cur_index, down_index) = - coe;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
            else
                B(cur_index, 1) = B(cur_index, 1) + heat_from_chip;
            end
            % top neighbor
            if (z ~= (c - 1))
                top_index = cur_index +  a * b;
                coe = k / size * area;
                A(cur_index, top_index) = - coe;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
            else
                coe = h * area;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
            end
        end
    end
end

% find the solution
T = A \ B;

T_plot = zeros(a, b, c);
for x = 1 : a
    for y = 1 : b
        for z = 1 : c
            index = (z - 1) * (a * b) + (y - 1) * a + x;
            T_plot(x, y, z) = T(index);
        end
    end
end

[x, y, z] = meshgrid(1 : a, 1 : b, 1 : c);
zslice = 1 : c;
slice(x, y, z, T_plot, [], [], zslice);
colormap hot;
colorbar;
% v = rand([10 10 10]);
% [x,y,z] = meshgrid(1:10,1:10,1:10);
% zslice = 1:10;
% slice(x,y,z,v,[],[],zslice)
% colormap hot;
% colorbar;