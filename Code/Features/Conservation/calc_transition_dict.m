function [trans_dict, probs] = calc_transition_dict(RNA)
%Helper function for prob, calculating transition dictionary and matrix based on RNA
%(Markov order 1)

%RNA = The entire mRNA strand

alphabet = {'a','c','g','u'};

keys = {};
index = 1;
for i = 1:length(alphabet)
    for j = 1:length(alphabet)
        keys{index} = [alphabet{i},alphabet{j}];
        index = index + 1;
    end
end

all_2mers = [regexp(RNA,'\w{2}','match') regexp(RNA(2:end),'\w{2}','match')];
all_2mers = sort(all_2mers);

probs = zeros(length(alphabet));

for i = 1:length(alphabet)
    for j = 1:length(alphabet)
        pair = [alphabet{i},alphabet{j}];
            probs(i,j) = sum(count(all_2mers,pair));
    end
    if sum(probs(i,:)) ~= 0
    probs(i,:) = probs(i,:) / sum(probs(i,:));
    end
end

trans_dict = containers.Map(keys,reshape(probs',1,length(keys)));

end