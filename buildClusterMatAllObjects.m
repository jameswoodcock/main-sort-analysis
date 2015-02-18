close all
clear all

pkg load io

experience = csvread('./data/experience.csv');
material = 'SS';
col = ['B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z'];
ind = 0;
countGroups = 0;
sortDataAllExp = [];
sortDataAllNonExp = [];
[allData,y,z] = odsread('20150218_all_objects_data_missing_vals.ods')

for n = 1:21
	if n == 5
		continue
	end
	if exist('groupLabs')
		clear groupLabs
	end
	n
	[x,y,z] = odsread(['./data/meta_labels_P',num2str(n),'.ods']);
	if ~isempty(x)
		sortData = zeros(length(x),length(allData));
		
		if experience(n+1,3) == 1
			for m = 1:max(x(:,1))
				groupLabs{m} = [char(y(m,1)),'_P',num2str(n)];	
			end

			for p = 1:length(allData)
				if allData(p,n) ~= 0
				sortData(allData(p,n),p) = 1;
				end
			end

			if exist('groupLabsAllExp')
				groupLabsAllExp = [groupLabsAllExp groupLabs];
			else
				groupLabsAllExp = groupLabs;
			end

			sortDataAllExp = [sortDataAllExp; sortData];
		else
			for m = 1:max(x(:,1))
				groupLabs{m} = [char(y(m,1)),'_P',num2str(n)];	
			end

			for p = 1:length(allData)
				if allData(p,n) ~= 0
				sortData(allData(p,n),p) = 1;
				end
			end

			if exist('groupLabsAllNonExp')
				groupLabsAllNonExp = [groupLabsAllNonExp groupLabs];
			else
				groupLabsAllNonExp = groupLabs;
			end

			sortDataAllNonExp = [sortDataAllNonExp; sortData];
		end

	end
end	%end n

groupLabsAllExp = groupLabsAllExp';
groupLabsAllNonExp = groupLabsAllNonExp';

save(['All.mat'], '-mat', 'sortDataAllExp','groupLabsAllExp', 'sortDataAllNonExp','groupLabsAllNonExp')



