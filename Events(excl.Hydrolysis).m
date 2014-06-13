function [s, b, l] = event(s_in, b_in)% s- матрица с димерами 13х5000 (2 - GTP, 1 - GDP, 0 - нет), b - матрица с латеральными связями 
(1 - есть, 0 - нет. Если это 13 протофиламент, то у него 2 - есть обе мономерные связи, 1 - одна(нижняя) мономерная связь, 0 - нет мономерных
связей), l - длина микротрубочки
%#codegen

minus_1 = zeros(1, 65000); %minus - отрыв димера. данный блок рассчитывает время для отрыва всех возможных димеров, и записывает это время 
и индекс димера
minus_2 = zeros(1, 65000);
minus_3 = zeros(1, 65000);
for idx = 1:numel(s_in),
    k = 0; %фиктивная переменная, обозначает, что есть хотя бы одна латеральная связь. Далее мы это проверяем, и если латеральных связей нет,
    заменяем k = 1
    if and((s_in(idx) > 0), (b_in(idx) == 0))%если есть димер на позиции, и нет латеральных связей справа
        if (idx - size(b_in, 2) > 0)%если это ПФ №2-13
            if (b_in(idx - size(b_in, 2)) == 0) %если лат. связь слева отсутствует 
                k = 1;
            end
        else% если это ПФ №1
            [row, col] = ind2sub(size(b_in), idx);% переводим порядковый номер элемента матрицы из числа в формат (строка, столбец)
            if (row > 2)%если этот димер не входит в первые два ряда внизу - там дырка, нет для него замыкающего димера на 13м ПФ
                if (and(b_in(row-2, 13) < 2, b_in(row-1, 13) == 0)); %нет обеих латеральных связей слева
                    k = 1;
                end
            end
        end
    end
    if k == 1% если нет обеих лат. связей 
        minus_1(find(minus_1==0, 1, 'first')) = 1;
        minus_2(find(minus_2==0, 1, 'first')) = idx;
        if (s_in(idx - 1) == 2)%если он соединен с GTP-Tu
            minus_3(find(minus_3==0, 1, 'first')) = (-log(rand))/0.02; % вероятность его отрыва
        else% если он соединен с GDP-Tu
            minus_3(find(minus_3==0, 1, 'first')) = (-log(rand))/20.00;% вероятность его отрыва
        end
    end
end

        
plus_1 = zeros(1, 65000); %вероятности для присоединения димеров
plus_2 = zeros(1, 65000);
plus_3 = zeros(1, 65000);
for idx = 1:numel(s_in)-1,
    if and((s_in(idx) > 0), (s_in(idx + 1) == 0)) %если есть димер на позиции и следом за ним - нет димера, к нему можно присоединить димер
        plus_1(find(plus_1 == 0, 1, 'first')) = 2;
        plus_2(find(plus_2 == 0, 1, 'first')) =  idx + 1;
        plus_3(find(plus_3 == 0, 1, 'first')) =  (-log(rand))/1.25;%вероятность присоединения
    end
end


break_lateral_1 = zeros(1, 65000); %вероятности разрыва латеральной связи из всех возможных
break_lateral_2 = zeros(1, 65000);
break_lateral_3 = zeros(1, 65000);
for idx = 1:numel(b_in),
    [row, col] = ind2sub(size(b_in), idx);
    if (and(b_in(idx) > 0, sum(b_in(row+1:size(b_in, 2), col)) == 0))%если выше нет лат. связей
        if (col == 13)%Рассматриваем 13й пф
            if (b_in(idx) == 1)%если только нижняя связь
                if (and(b_in(row+1, 1) == 1, b_in(row, 12) == 1))%если есть связи справа и слева
                    pi = 1000;
                else
                    pi = 1;
                end
            else%есть обе связи
                if (and(b_in(row+2, 1) == 1, b_in(row, 12) == 1))%если есть связи справа и слева
                    pi = 1000;
                else
                    pi = 1;
                end
            end
            if (and(s_in(idx) == 1, s_in(row+1,1) == 1))%если оба димера, образующих эту связь - GDP
                k_break = 800;
            elseif (or(and(s_in(idx) == 1, s_in(row+1,1) == 2), and(s_in(idx) == 2, s_in(row+1,1) == 1)))%если один GDP, а второй GTP
                k_break = 180;
            else%оба GTP
                k_break = 140;    
            end    
        else
            if (col == 1) %Первый ПФ
                if (and(and(b_in(row-2, 13) == 2, b_in(row, 2) == 1), b_in(row-1, 13) >= 1))%если вокруг данной связи есть все лат. связи
                % две - на 13м ПФ и одна - на 2м пф, уменьшаем пи-брейк в 1000 раз
                    pi = 1000;
                else
                    pi = 1;
                end
            else
                if(and(b_in(row, col-1) == 1, b_in(row, col+1) == 1))% для остальных ПФ лат. связи окружающих - две. Справа и слева. Если они есть, уменьшаем пи-брейк в 1000 раз
                    pi = 1000;
                else
                    pi = 1;
                end
            end
            if (and(s_in(idx) == 1, s_in(row,col+1) == 1))%если оба GDP
                k_break = 400;
            elseif (or(and(s_in(idx) == 1, s_in(row,col+1) == 2), and(s_in(idx) == 2, s_in(row,col+1) == 1)))%если один GDP, а второй GTP
                k_break = 90;
            else% both GTP
                k_break = 70;
            end
        end
    break_lateral_1(find(break_lateral_1 == 0, 1, 'first')) = 3;
    break_lateral_2(find(break_lateral_2 == 0, 1, 'first')) = idx;
    break_lateral_3(find(break_lateral_3 == 0, 1, 'first')) = (-log(rand))/(k_break/pi);
    end
end
            
        
bond_lateral_1 = zeros(1, 65000); %вероятности образования лат. связей
bond_lateral_2 = zeros(1, 65000);
bond_lateral_3 = zeros(1, 65000);
for idx = 1:numel(b_in),
    [row, col] = ind2sub(size(b_in), idx);
    if col < 13  % ПФ №1-12
        if (and(and(b_in(idx) == 0, s_in(idx) > 0), and(s_in(idx+size(b_in, 2)) > 0, sum(b_in(1:row-1, col) == 0) == 0))) %если внизу есть все лат связи
            bond_lateral_1(find(bond_lateral_1 == 0, 1, 'first')) = 4;
            bond_lateral_2(find(bond_lateral_2 == 0, 1, 'first')) = idx;
            bond_lateral_3(find(bond_lateral_3 == 0, 1, 'first')) =  (-log(rand))/100;
        end
    else%13 ПФ
        if (and(b_in(idx) < 2, sum(b_in(1:row-1, col) == 0) == 0)),%если внизу есть лат связи
            if or(and(and((b_in(idx) == 0), s_in(idx) > 0), s_in(row-1, 1) > 0), (and(and((b_in(idx) == 1), s_in(idx) > 0), s_in(row-2, 1) > 0)))%если есть димеры для предполагаемой связи, а самой связи нет
                bond_lateral_1(find(bond_lateral_1 == 0, 1, 'first')) = 4;
                bond_lateral_2(find(bond_lateral_2 == 0, 1, 'first')) = idx;
                bond_lateral_3(find(bond_lateral_3 == 0, 1, 'first')) =  (-log(rand))/100;%есть только нижняя мономерная связь
            end
        end
    end
end

final_list_1 = [minus_1, plus_1, break_lateral_1, bond_lateral_1]; %находим меньшее время исполнения
final_list_2 = [minus_2, plus_2, break_lateral_2, bond_lateral_2];
final_list_3 = [minus_3, plus_3, break_lateral_3, bond_lateral_3];
final_list_3(final_list_3 == 0) = inf;
least_time = find(final_list_3 == min(final_list_3));
least_time_1 = least_time(randi(numel(least_time)));
idx = final_list_2(least_time_1);
if (final_list_1(least_time_1) == 1)  % выполняем данное событие
    s_in(idx) = 0;
elseif (final_list_1(least_time_1) == 2)
    s_in(idx) = 2;
elseif (final_list_1(least_time_1) == 3)
    b_in(idx) = b_in(idx) -  1;
else
    b_in(idx) = b_in(idx) +  1;
end


l = nnz(s_in)/13 - 10;     %считаем длину МТ 
s = s_in;
b = b_in;
