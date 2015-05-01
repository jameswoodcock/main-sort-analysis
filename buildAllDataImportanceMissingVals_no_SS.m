close all
clear all

pkg load io


sheet = ['RD';'FF';'ND';'LE'];
col = ['B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z'];
ind = 0;
for n = 1:21
objects = [];
for m = 1:4
n
[x,y,z] = odsread(['./data/all_objects_P',num2str(n),'.ods'],sheet(m,:));
if isempty(x)
objects = [objects;nan(length(y),1)];
else
objects = [objects;x(:,3)];
end

if m == 4
ind = ind + 1;
odswrite('all_objects_importance_missing_vals.ods',objects,'Sheet1',[col(n),num2str(1),':',col(n),num2str(length(objects))]);
save all_objects_importance_missing_vals_no_SS.mat
end

end
end
