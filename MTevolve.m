function MTevolve(N,M)
l=zeros(1,M);
t=zeros(1,M);
s_in = zeros(N, 13);
b_in = zeros(N, 13);


[s_in, b_in, l(1),t(1), sob(1), idx_final(1)]=event(N, s_in, b_in);

filename = '1.gif';
for i=2:M
      [s_in, b_in, l(i),dt, sob(i), idx_final(i)]=event(N, s_in, b_in);
      t(i)=t(i-1)+dt;
      if (rem(i, 1000) == 0)
          imagesc(s_in);
%           textStrings = num2str(b_in(:),'%d');
%           textStrings = strtrim(cellstr(textStrings));
%           [x,y] = meshgrid(1:13, 1:N);
%           hStrings = text(x(:),y(:),textStrings(:),...      
%                     'HorizontalAlignment','center');
%           set(hStrings,'Color','white');
          frame = getframe();
          im = frame2im(frame);
          [imind,cm] = rgb2ind(im,256);
          if i == 1000;
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
          else
              imwrite(imind,cm,filename,'gif','WriteMode','append');
          end
      end
end

figure
plot(t,l,'r.-');
% for i=1:M
%     if (and(idx_final(i) > 6000, s_in(idx_final(i)) > 0))
%         if s_in(idx_final(i)+1) == 0
%             plot(i, sob(i), 'b.');
%             hold on
%         end
%     end
% end
