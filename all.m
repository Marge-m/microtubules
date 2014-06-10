%создаю две матрицы - одна с димерами(2 - GTP, 1 - GDP, 0 - пусто), одна со связями(имеется в виду связь с правым соседом). для ПФ 1-12 связь может быть 
1(есть) или 0(нет), а для 13-го 2(обе мономерных), 1(одна мономерная которая ближе к основанию МТ), 0(нет связей)


function [s, b] = mat_init(s_in, b_in)
%#codegen

if sum(s_in(:)) == 0 
  s_in(1:10,:) = 2; %GTP - Tu - затравка
end

if sum(b_in(:)) == 0
  b_in(1:10, 1:12) = 1; % димерные связи
  b_in(1:10, 13) = 2; % мономерные связи 
% для 13-го ПФ: 0 - нет связей, 1 - есть одна нижняя связь, 2 - есть обе связи
end

s = s_in;
b = b_in;





%описываю возможные события


function [s, b, l] = event(s_in, b_in)%k - длина МТ
%#codegen

minus = []; %события для отрыва димера
for idx = 1:numel(s_in),
    k = 0; фиктивная переменная
    if and((s_in(idx) > 0), (b_in(idx) == 0))%если есть димер на позиции, и нет латеральных связей справа
        if (idx - size(b_in, 2) > 0)%если это ПФ №2-13
            if (b_in(idx - size(b_in, 2)) == 0) %если лат. связь слева отсутствует 
                k = 1;
            end
        else% если это ПФ №1
            [row, col] = ind2sub(size(b_in), idx);
            if (row > 2)%если этот димер не входит в первые два ряда внизу - там дырка, нет для него замыкающего димера на 13м ПФ
                if (and(b_in(row-2, 13) < 2, b_in(row-1, 13) == 0)); %нет обеих мономерных латеральных связей слева
                    k = 1;
                end
            end
        end
    end
    if k == 1% если нет обеих лат. связей 
        if (s_in(idx - 1) == 2)%если он соединен с GTP-Tu
            minus = [minus; ['minus' idx (-log(rand))/0.02]];% вероятность его отрыва
        else% если он соединен с GDP-Tu
            minus = [minus; ['minus' idx (-log(rand))/20.00]];% вероятность его отрыва
        end
    end
end

        
plus = [] %присоединение димера
for idx = 1:numel(s_in),
    if and((s_in(idx) > 0), (s_in(idx + 1) == 0)) %если есть димер на позиции и следом за ним - нет димера
        plus = [plus; ['plus' (idx + 1) (-log(rand))/1.25]];%вероятность присоединения
    end
end


break_lateral = [] %разрыв латеральной связи
for idx = 1:numel(b_in),
    [row, col] = ind2sub(size(b_in), idx)
    if (and(b_in(idx) > 0, sum(b_in(row+1:size(b_in, 2), col)) == 0))%если выше нет лат. связей
        if (col == 13)%13 пф
            if (b_in(idx) == 1)%если есть только нижняя мономерная связь у данного димера, рассматриваем её - действует ли пи-брейк
                if (and(b_in(row+1, 1) == 1, b_in(row, 12) == 1))%если есть связи справа и слева от данной латеральной связи
                    pi = 1000;
                else %если связей нет
                    pi = 1;
                end
            else% если есть обе мономерных связи, рассматриваем верхнюю - действует ли пи-брейк
                if (and(b_in(row+2, 1) == 1, b_in(row, 12) == 1))%если есть связи справа и слева
                    pi = 1000;
                else
                    pi = 1;
                end
            end
            if (and(s_in(idx) == 1, s_in(row+1,1) == 1))%если оба GDP (отрывающийся и тот, от которого отрывается)
                k_break = 800;
            elseif (or(and(s_in(idx) == 1, s_in(row+1,1) == 2), and(s_in(idx) == 2, s_in(row+1,1) == 1)))%если один GDP, а второй GTP
                k_break = 180;
            else %оба GTP
                k_break = 140;    
            end    
        else
            if (col == 1) % для первого ПФ
                if (and(b_in(row-2, 13) == 2, b_in(row, 2) == 1, b_in(row-1, 13) >= 1)) %смотрим есть ли латеральные связи слева дальной связи
                    pi = 1000;
                else % если их нет
                    pi = 1;
                end
            else % для ПФ № 2-12
                if(and(b_in(row, col-1) == 1, b_in(row, col+1) == 1)) %если рядом есть лат. связт
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
    break_lateral = [break_lateral; ['break_lat' (idx) (-log(rand)/(k_break/pi))]];
    end
end
            
        
bond_lateral = []; %образование латеральной связи справа данной позиции
for idx = 1:numel(b_in),
    [row, col] = ind2sub(size(b_in), idx)
    if col < 13   % ПФ №1-12
        if (and(b_(idx) == 0, s_in(idx) > 0, s_in(idx+size(b_i, 2)) > 0, sum(b_in(1:row-1, col) == 0) == 0)) %если внизу есть все лат связи, образование связи возможно
            bond_lateral = [bond_lateral; [('bond_lat') (idx) (col) (-log(rand)/100)]];
        end
    else %13 ПФ
        if (and(b_in(idx) < 2, sum(b_in(1:row-1, col) == 0) == 0)),%если внизу есть лат связи, образование связи возможно
            if (and((b_in(idx) == 0), s_in(idx) > 0, s_in(row-1, 1) > 0))%если есть димеры для предполагаемой связи, а самой мономерной связи нет (ни нижней, ни верхней)
                bond_lateral = [bond_lateral; [('bond_lat') (idx) (-log(rand)/100)]];%есть только нижняя мономерная связь
            end
            if (and((b_in(idx) == 1), s_in(idx) > 0, s_in(row-2, 1) > 0))%если есть димеры для предполагаемой связи, а самой связи нет (верхней)
                bond_lateral = [bond_lateral; [('bond_lat') (idx) (-log(rand)/100)]];%есть обе мономерные связи
            end
        end
    end
end

final_list = [minus; plus; break_lateral; bond_lateral]; %соединили в матрицу события
least_time = find(final_list(:,3) == min(final_list(:,3))); % рассчитали меньшее время 
idx = final_least(2, least_time) % нашли индекс димера, над которыммм совершается событие
[row, col] = ind2sub(size(b_in), idx);
if (final_list(1, least_time) == 'minus') % убираем димер
    s_in(idx) = 0;
elseif (final_list(1, least_time) == 'plus') % присоединяем димер
    s_in(idx) = 2;
elseif (final_list(1, least_time) == 'break_lat') % разрываем лат. связь
    b_in(idx) = b_in(idx) -  1;
else % создаем лат. связь
    b_in(idx) = b_in(idx) +  1;
end

l = nnz(s_in)/13 - 10; %считаем длину МТ минус 10 (затравка). После чего подаем эту длину для построения графика

        
s = s_in;
b = b_in;


% Совершаем гидролиз

function s = hydr(s_in)
%#codegen

ind = find([zeros(10, size(s_in, 2)); s_in(11:5000, :)] == 2); % находим все GTP-Tu кроме затравки 
[sel, c] = max( s_in == 0, [], 1); % находим первый нулевой элемент в каждом столбце
c1 = (c - 1) + ((0:12) * size(s_in, 1)); %  находим последний димер в каждом столбце
for item = c1,
    ind = ind(ind~=item); %убираем из списка  последние димеры
end 
 
s_in(ind(randi(numel(ind)))) = 1; % меняем рандомный ГТФ на ГДФ
s = s_in;
end


%цикл повторяется снова и на выход подается длина МТ, строится график
