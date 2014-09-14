function [s, b, l,t, sob, idx_final] = event(N, s_in, b_in)%  s- matrix with dimers(2 - GTP, 1 - GDP, 0 - none), b - matrix with lateral bonds (bond with a right neighbore) 
%(1 - +, 0 - none. For protofilament 13: 2 - both monomer bonds, 1 - only lower monomer bond, 0 - none
% l - length of MT
%#codegen
% input
% matrices-------------------------------------------------------------------




if sum(s_in(:)) == 0 
  s_in(1:10,:) = 2; %GTP - Tu - zatravka
end

if sum(b_in(:)) == 0
  b_in(1:10,1:12) = 1;% bonds for 1-12 PF
  b_in(1:10, 13) = 2;%bonds for 13th PF
end
%-------------------------------------------------------------------

k_minus_GTP = 0.02;
k_minus_GDP = 20;
k_plus = 12.5;
k_break1 = 800;
k_break2 = 180;
k_break3 = 140;
k_break4 = 400;
k_break5 = 90;
k_break6 = 70;
k_bond = 100;
k_hydr = 0.70;

minus_1 = zeros(1, N*13); %array for event name(code). '1' is for shorten
minus_2 = zeros(1, N*13);%2 - index of dimer
minus_3 = zeros(1, N*13);%3 -execution time
for idx = 1:numel(s_in),
    k = 0; %dummy variable
    if and((s_in(idx) > 0), (b_in(idx) == 0))%if dimer exists and there is no right lat. bond
        if (idx - size(b_in, 1) > 0)%PF 2-13
            if (b_in(idx - size(b_in, 1)) == 0) %if there is no left lat. bond
                k = 1;
            end
        else% PF 1
            if (rem(idx-1, size(b_in, 1))+1 > 2)%row > 2 (this is evaluation of row. Index is given, we transform it to view (row, column)
                if (and(b_in(rem(idx-1, size(b_in, 1))+1-2, 13) < 2, b_in(rem(idx-1, size(b_in, 1))+1-1, 13) == 0)); %if thre is no monomer bonds to the left
                    k = 1;
                end
            end
        end
    end
    if k == 1% if dimer is free from lat. bonds
        minus_1(find(minus_1==0, 1, 'first')) = 1;% event 
        minus_2(find(minus_2==0, 1, 'first')) = idx;
        if (s_in(idx - 1) == 2)%if it is bound to GTP-Tu
            minus_3(find(minus_3==0, 1, 'first')) = (-log(rand))/k_minus_GTP; % exec. time 0.02
        else% bound to GDP-Tu
            minus_3(find(minus_3==0, 1, 'first')) = (-log(rand))/k_minus_GDP;% exec. time 20
        end
    end
end

        
plus_1 = zeros(1, N*13); % k+
plus_2 = zeros(1, N*13);
plus_3 = zeros(1, N*13);
for idx = 1:numel(s_in)-1,
    if and((s_in(idx) > 0), (s_in(idx + 1) == 0)) %if dimer exists and is last in PF
        plus_1(find(plus_1 == 0, 1, 'first')) = 2;
        plus_2(find(plus_2 == 0, 1, 'first')) =  idx + 1;
        plus_3(find(plus_3 == 0, 1, 'first')) =  (-log(rand))/k_plus;%exec time 1.25
    end
end


break_lateral_1 = zeros(1, N*13); %break of lat. bond
break_lateral_2 = zeros(1, N*13);
break_lateral_3 = zeros(1, N*13);


for idx = 1:numel(b_in),
    row = rem(idx-1, size(b_in, 1))+1;
    col = (idx-row)/size(b_in, 1) + 1;
    if (and(row<size(b_in,1), b_in(idx) > 0))%we don't take into account the last row. And check if dimer has lat. bond
        if (sum(b_in(row+1:size(b_in, 1), col)) == 0)%if all lat. bonds below exists
            if (col == 13)%PF 13
                if (b_in(idx) == 1)%if dimer has only one monomer bond
                    if (and(b_in(row+1, 1) == 1, b_in(row, 12) == 1))%if neighbore lat. bonds exists, we decrease K by pi
                        pi = 1000;
                    else
                        pi = 1;
                    end
                else%if both monomer bonds
                    if (and(b_in(row+2, 1) == 1, b_in(row, 12) == 1))%neighbore bonds exists
                        pi = 1000;
                    else
                        pi = 1;
                    end
                end
                if (and(s_in(idx) == 1, s_in(row+1,1) == 1))%if both dimers which are making a bond are GDP
                    k_break = k_break1;%800
                elseif (or(and(s_in(idx) == 1, s_in(row+1,1) == 2), and(s_in(idx) == 2, s_in(row+1,1) == 1)))%one is GDP, the other is GTP
                    k_break = k_break2;%180
                else%both GTP
                    k_break = k_break3; %140   
                end    
            elseif (col == 1)
                if row > 2%PF 1
                    if (and(and(b_in(row-2, 13) == 2, b_in(row, 2) == 1), b_in(row-1, 13) >= 1))%neighbore lat. bonds are present
                        pi = 1000;
                    else
                        pi = 1;
                    end
                    if (and(s_in(idx) == 1, s_in(row,col+1) == 1))%both GDP
                        k_break = k_break4; %400
                    elseif (or(and(s_in(idx) == 1, s_in(row,col+1) == 2), and(s_in(idx) == 2, s_in(row,col+1) == 1)))%one GDP, the other GTP
                        k_break = k_break5;%90
                    else% both GTP
                        k_break = k_break6;%70
                    end
                end
            elseif (col == 12)%PF 12
                if (and(b_in(row, 13) == 2, b_in(row, 11) == 1))%neighbore lat. bonds are present
                    pi = 1000;
                else
                    pi = 1;
                end
                if (and(s_in(idx) == 1, s_in(row, 13) == 1))%both GDP
                    k_break = k_break4;%400
                elseif (or(and(s_in(idx) == 1, s_in(row, 13) == 2), and(s_in(idx) == 2, s_in(row, 13) == 1)))%one GDP, the other GTP
                    k_break = k_break5;%90
                else% both GTP
                    k_break = k_break6;%70
                end
            else%PF 2-11
                if(and(b_in(row, col-1) == 1, b_in(row, col+1) == 1))% both neighbore lateral bonds are present
                    pi = 1000;
                else
                    pi = 1;
                end
                if (and(s_in(idx) == 1, s_in(row,col+1) == 1))%both GDP
                    k_break = k_break4;%400
                elseif (or(and(s_in(idx) == 1, s_in(row,col+1) == 2), and(s_in(idx) == 2, s_in(row,col+1) == 1)))%one GDP, the other GTP
                    k_break = k_break5;%90
                else% both GTP
                    k_break = k_break6;%70
                end
            end
            if row > 2
                break_lateral_1(find(break_lateral_1 == 0, 1, 'first')) = 3;
                break_lateral_2(find(break_lateral_2 == 0, 1, 'first')) = idx;
                break_lateral_3(find(break_lateral_3 == 0, 1, 'first')) = (-log(rand))/(k_break/pi);%exec. time
            end
        end
    end
end        
 


bond_lateral_1 = zeros(1, N*13); %bond formation
bond_lateral_2 = zeros(1, N*13);
bond_lateral_3 = zeros(1, N*13);
for idx = 1:numel(b_in),
    row = rem(idx-1, size(b_in, 1))+1;
    col = (idx-row)/size(b_in, 1) + 1;
    if col < 13  % ПФ 1-12
        if (and(b_in(idx) == 0, s_in(idx) > 0))%we have a dimer but it is not laterally bound to a right neighbore
            if (and(s_in(idx+size(b_in, 1)) > 0, sum(b_in(1:row-1, col) == 0) == 0)) %right neighbore is present and between these two PF all below lat. bonds are present
                bond_lateral_1(find(bond_lateral_1 == 0, 1, 'first')) = 4;
                bond_lateral_2(find(bond_lateral_2 == 0, 1, 'first')) = idx;
                bond_lateral_3(find(bond_lateral_3 == 0, 1, 'first')) =  (-log(rand))/k_bond;%exec. time 100
            end
        end
    else%13 PF
        if (b_in(idx) < 2)%less than both monomer bonds are present
            if (sum(b_in(1:row-1, col) == 0) == 0)%below all lat. bonds are present
                if row > 2
                    if or(and(and((b_in(idx) == 0), s_in(idx) > 0), s_in(row-1, 1) > 0), (and(and((b_in(idx) == 1), s_in(idx) > 0), s_in(row-2, 1) > 0)))%first case - no monomer bonds, second - one monomer bond
                        bond_lateral_1(find(bond_lateral_1 == 0, 1, 'first')) = 4;
                        bond_lateral_2(find(bond_lateral_2 == 0, 1, 'first')) = idx;
                        bond_lateral_3(find(bond_lateral_3 == 0, 1, 'first')) =  (-log(rand))/k_bond;%exec time
                    end
                end
            end
        end
    end
end

hydr_1 = zeros(1, N*13); %hydrolysis
hydr_2 = zeros(1, N*13);
hydr_3 = zeros(1, N*13);
for idx = 1:numel(s_in)
    row = rem(idx-1, size(b_in, 1))+1;%row
    if row > 10%excluding zatravka
        if s_in(idx) == 2%dimer is GTP-Tu
            if (s_in(idx+1) ~= 0)%excluding last row of GTP-Tu dimers
                hydr_1(find(hydr_1 == 0, 1, 'first')) = 5;
                hydr_2(find(hydr_2 == 0, 1, 'first')) = idx;
                hydr_3(find(hydr_3 == 0, 1, 'first')) = (-log(rand))/k_hydr;%exec.time
            end
        end
    end
end


final_list_1 = [minus_1, plus_1,  bond_lateral_1, break_lateral_1, hydr_1]; %list of events
final_list_2 = [minus_2, plus_2, bond_lateral_2, break_lateral_2, hydr_2];%list of indexes
final_list_3 = [minus_3, plus_3,  bond_lateral_3, break_lateral_3, hydr_3];%list of exec. time
final_list_3(final_list_3 == 0) = inf;%exclude empty spaces in array of exec. time 
least_time = find(final_list_3 == min(final_list_3));%find minimal time
least_time_1 = least_time(randi(numel(least_time)));%find event of minimal time
idx_final = final_list_2(least_time_1);%index
sob1 = final_list_1(final_list_1 ~= 0);
sob2 = final_list_2(final_list_2 ~= 0);
sob3 = final_list_3(final_list_3 ~= inf);
sob = final_list_1(least_time_1);

if (sob == 1)  % first event - shorten
    s_in(idx_final:(ceil(idx_final/size(s_in, 1)) * size(s_in, 1))) = 0;
elseif (sob == 2)%etc.
    s_in(idx_final) = 2;
elseif (sob == 3)
    b_in(idx_final) = b_in(idx_final) -  1;
elseif (sob == 4)
    b_in(idx_final) = b_in(idx_final) +  1;
else
    s_in(idx_final) = 1;
end

t=final_list_3(least_time_1);%least time to plot
l = nnz(s_in)/13 - 10;     % length of MT 
s = s_in;
b = b_in;
