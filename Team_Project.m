%% ME 3345
%% Spring 2015
%% Term project -- Topic 2 (Heat sink problem)
%% AUTHOR:  Tiancheng Gong (tgong7)

%% PURPOSE: 
%      Find and plot the temperature distribution of the heat sink, find
%  the position with the highest temperature using numerical method.

%% ASSUMPTIONS: 
%      - The heat transfer coefficient is uniform.
%      - The heat flux provided by the generation in the chip is uniform.
%      - Steady state.
%      - 3D conduction.
%      - No radiation happens.
%      - The size of the block unit used is .5mm*.5mm*.5mm. The temperature
%      within in each block is uniform.
%      - The temperature of the heat sink is around 300K.
%      - The temperature of the infinite distance is around 300K. It is the
%      difference between the environment and the chip instead of the
%      temperature of environment itself that matters in this problem.

%% IMPORTANT VALUES:
%      k = 237 W/m*K;
%      h = 100 W/m^2*K;
%      size = 0.5 mm = 5e(-4) m;
%      Total power = 4W;

%% ANALYSIS AND METHOD:
%      - Since the shape of the heat sink is symmetric, the model can be
%        simplified  by dividing the heat sink into 4 identical parts, 4
%        squares and one rod in the middle of each square. (It can be
%        divided into at most 8 identical triangle parts, but it is not
%        easy to build this model using cube blocks.)
%      - The tool used here is the finite difference equations. For each
%        cube block unit, we can derive a finite difference equation based
%        on the temperatrue of its neighbors and the property of the
%        material. We will have n cubes as n unknowns as well as n
%        equations, thus the equations are solvable.
%      - Using the solution derived, we can thus plot the temperature 
%        distribution and find the position of the heighest temperature.


%% PART 1: Create the array for the coefficients of the equations.
%  Linear Equations: At = B; B is all zeros.
size = .0005;
a = .008 / size;
b = .008 / size;
c = .006 / size;
l = .02 / size;
blocks_in_base = a * b * c;
blocks_in_rod = 32 * l;
total_blocks = blocks_in_base + blocks_in_rod;
A = zeros(total_blocks, total_blocks); % For all blocks
B = zeros(total_blocks, 1);

% Set the value of coefficients based on finite difference theorem.
% Consider the border cubes carefully.
area = size ^ 2;
k = 237;
h = 100;
total_power = 4;
rm_temp = 300;
heat_from_chip = total_power / (4 *(a * b));
% the top plane for the base 0: nothing above it; 0.5: triangle rod block
% is above it; 1: cube rod block is above it
top_plane = zeros(a, b);
top_plane(6:11, 7:10) = .5;
top_plane(7:10, 6:11) = .5;
top_plane(6:11, 8:9) = 1;
top_plane(7:10, 7:10) = 1;
top_plane(8:9, 6:11) = 1;

% collect all the finite difference equations of the bese part.
count = 0; % count the top blocks connecting to the rod blocks
for z = 0 : (c - 1)
    for y = 0 : (b - 1)
        for x = 0 : (a - 1)
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
                B(cur_index, 1) = B(cur_index, 1) + coe * rm_temp;
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
                B(cur_index, 1) = B(cur_index, 1) + coe * rm_temp;
            end
            % back neighbor
            if (y ~= (b - 1))
                back_index = cur_index + a;
                coe = k / size * area;
                A(cur_index, back_index) = - coe;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
            end
            % bot neighbor
            if (z ~= 0)
                bot_index = cur_index - a * b;
                coe = k / size * area;
                A(cur_index, bot_index) = - coe;
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
                % check whether the block is connected to the rod
                % top neighbor is the environment
                if (top_plane(x+1, y+1) == 0)
                    coe = h * area;
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                    B(cur_index, 1) = B(cur_index, 1) + coe * rm_temp;
                % top neighbor is the cube rod block
                elseif (top_plane(x+1, y+1) == 1)
                    count = count + 1;
                    top_index = blocks_in_base + count;
                    coe = k / size * area;
                    A(cur_index, top_index) = - coe;
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                % top neighbor is the triangle rod block
                else
                    count = count + 1;
                    top_index = blocks_in_base + count;
                    coe_k = k / size * area / 2;
                    coe_h = h * area / 2;
                    A(cur_index, top_index) = - coe_k;
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe_k + coe_h;
                    B(cur_index, 1) = B(cur_index, 1) + coe_h * rm_temp;
                end
            end
        end
    end
end


count = 0;
index_matrix = zeros(6, 6);
condition = zeros(6, 6);
for j = 1 : 6
    for i = 1: 6
        if (top_plane(i + 5, j + 5) ~= 0)
            condition(i, j) = top_plane(i + 5, j + 5);
            count = count + 1;
            index_matrix(i, j) = count;
        end
    end
end
side1 = [2, 3, 11, 16, 17, 22, 30, 31];
side2 = [1, 4, 5, 10, 23, 28, 29, 32];

for z = 1 : l
    for y = 1 : 6
        for x = 1 : 6
            if (condition(x, y) ~= 0)
                cur_index = blocks_in_base + (z - 1) * 32 + index_matrix(x, y);
                % left neighbor
                if ((x ~= 1) && (condition(x - 1, y) ~= 0))
                    left_index = cur_index - 1;
                    coe = k / size * area;
                    A(cur_index, left_index) = - coe;
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                end
                % right neighbor
                if ((x ~= 6) && (condition(x + 1, y) ~= 0))
                    right_index = cur_index + 1;
                    coe = k / size * area;
                    A(cur_index, right_index) = - coe;
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                end
                % front neighbor
                if ((y ~= 1) && (condition(x, y - 1) ~= 0))
                    front_index = blocks_in_base + (z - 1) * 32 + index_matrix(x, y - 1);
                    coe = k / size * area;
                    A(cur_index, front_index) = - coe;
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                end
                % back neighbor
                if ((y ~= 6) && (condition(x, y + 1) ~= 0))
                    back_index = blocks_in_base + (z - 1) * 32 + index_matrix(x, y + 1);
                    coe = k / size * area;
                    A(cur_index, back_index) = - coe;
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                end
                % top neighbor
                if (z ~= l)
                    top_index = cur_index + 32;
                    coe = k / size * area * condition(x, y); % whether the cross section is a triangle
                    A(cur_index, top_index) = - coe;
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                else
                    coe = h * area * condition(x, y);
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                    B(cur_index, 1) = B(cur_index, 1) + coe * rm_temp;
                end
                % bot neighbor
                if (z ~= 1)
                    bot_index = cur_index - 32;
                else
                    bot_index = blocks_in_base - a * b + (y + 4) * a + x + 5;
                end
                coe = k * condition(x, y)/ size * area;
                A(cur_index, bot_index) = - coe;
                A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                % consider the h
                if (~isempty(find(side1 == index_matrix(x, y), 1)))
                    coe = h * area;
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                    B(cur_index, 1) = B(cur_index, 1) + coe * rm_temp;
                elseif (~isempty(find(side2 == index_matrix(x, y), 1)))
                    coe = h * area * 1.41421;
                    A(cur_index, cur_index) = A(cur_index, cur_index) + coe;
                    B(cur_index, 1) = B(cur_index, 1) + coe * rm_temp;
                end
            end
        end
    end
end

%% PART 2: Find the solution
T = A \ B;


%% PART 3: Plot and interpret the solution.
% put the temperatre soulutions into T_plot matrix
% T_plot = ones(a, b, c + l) * rm_temp;
T_base_plot = zeros(a, b, c);
T_rod_plot = ones(a, b, l) * rm_temp;
% put in the base T elements
for z = 1 : c
    for y = 1 : b
        for x = 1 : a
            index = (z - 1) * (a * b) + (y - 1) * a + x;
            T_base_plot(x, y, z) = T(index);
        end
    end
end


% plug in the rod T elements
index = blocks_in_base;
for z = 1 : l
    for y = 1 : b
        for x = 1 : a
            if (top_plane(x, y) ~= 0)
                index = index + 1;
                T_rod_plot(x, y, z) = T(index);
            end
        end
    end
end

figure(1)
[x, y, z] = meshgrid(0.5 : a - 0.5, .5 : b - .5, .5 : c - .5);
xslice = .5 : 1 : a - .5;
slice(x, y, z, T_base_plot, xslice, [], []);
colorbar;

figure(2)
[x, y, z] = meshgrid(5.5:13.5, 5.5:13.5, 13:l+12);
xslice = 6:2:13;
zslice = 13: 2: l+12;
slice(x, y, z, T_rod_plot(5:13, 5:13, :), xslice, [], zslice);
axis equal
colorbar;


figure(3)
[x, y, z] = meshgrid(0.5 : a - 0.5, .5 : b - .5, .5 : c - .5);
xslice = .5 : 1 : a - .5;
slice(x, y, z, T_base_plot, xslice, [], []);
hold

[x, y, z] = meshgrid(4.5:12.5, 4.5:12.5, 12:l+11);
xslice = 5:2:12;
zslice = 12: 2: l+11;
slice(x, y, z, T_rod_plot(5:13, 5:13, :), xslice, [], zslice);
axis equal
colorbar;

[T_max, max_loc] = max(T);
[T_min, min_loc] = min(T);
T_chip = T_max + .001 / (.008 * .008);
