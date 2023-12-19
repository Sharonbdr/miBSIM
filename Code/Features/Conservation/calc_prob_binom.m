function prob_binom = calc_prob_binom(RNA,RNA_start,seed_length)
%Calculating the binomial probability of the target site

%RNA = The whole RNA strand
%RNA_start = Start position of the seed in the RNA
%seed_length = Length between 6 and 7

RNA = lower(RNA);
RNA(RNA == 't') = 'u';
target_site = RNA(RNA_start:RNA_start + seed_length - 1);
[trans_dict,~] = calc_transition_dict(RNA);

%Calculating the probability of seeing the target sites given trans_dict

target_site_2mers = [regexp(target_site,'\w{2}','match') regexp(target_site(2:end),'\w{2}','match')];
probs = cellfun(@(x) trans_dict(x), target_site_2mers);
p = prod(probs);

%Calculating the binomial probability of seeing the target k times
k = count(RNA,target_site);
n = length(RNA) - seed_length + 1;

prob_binom = 1 - binocdf(k,n,p);
end