close all
clear all

pkg load io


sheet = ['RD';'FF';'ND';'SS';'LE'];
col = ['B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N'];
ind = 0;
for n = 1:10
objects = [];
for m = 1:5
n
[x,y,z] = odsread(['./data/all_objects_P',num2str(n),'.ods'],sheet(m,:));
if isempty(x)
break
end
objects = [objects;x(:,2)];

if m == 5
ind = ind + 1;
odswrite('all_objects_full_list.ods',objects,'Sheet1',[col(ind),num2str(1),':',col(ind),num2str(624)]);
end

end
end
