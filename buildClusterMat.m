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

for n = 1:21
	if exist('groupLabs')
		clear groupLabs
	end
	n
	[x,y,z] = odsread(['./data/all_objects_P',num2str(n),'.ods'],material);
	if ~isempty(x)
		sortData = zeros(max(x(:,1)),length(x));
		
		if experience(n+1,3) == 1
			for m = 1:max(x(:,1))
				groupLabs{m} = [char(y(m,10)),'_P',num2str(n)];	
			end

			for p = 1:length(x)
				sortData(x(p,1),p) = 1;
			end

			if exist('groupLabsAllExp')
				groupLabsAllExp = [groupLabsAllExp groupLabs];
			else
				groupLabsAllExp = groupLabs;
			end

			sortDataAllExp = [sortDataAllExp; sortData];
		else
			for m = 1:max(x(:,1))
				groupLabs{m} = [char(y(m,10)),'_P',num2str(n)];	
			end

			for p = 1:length(x)
				sortData(x(p,1),p) = 1;
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

save([material,'.mat'], '-mat', 'sortDataAllExp','groupLabsAllExp', 'sortDataAllNonExp','groupLabsAllNonExp')



