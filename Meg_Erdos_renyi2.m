%% ER MODEL
% 
clc;
close all;

% this sets the random seed according to the system clock, runs are now
% random
rng("shuffle")

% Abrams-Strogatz parameters
a = 1.31;
s = 0.6;   % status set so that language 1 (marked 0) is most popular

% Erdos-Renyi parameters
N = 1000;
p = 0.7;

init_ppl = binornd(1,0.5,1,N);  % initial condition with roughly half of each language
pct_ic = sum(init_ppl)/N;  % percent of IC speaking language 2
sim_max = 10;  % total number of simulations
t_max = 30;  % length of simulations

l1 = zeros(t_max,sim_max);  % fraction of people marked 0 in language vector
l2 = zeros(t_max,sim_max);  % fraction of people marked 1 in language vector

tic
for sim = 1:sim_max
    [G,n,m] = create_ER_Graph(N,p,sim,1);  % must use format 1 to have full matrix
    ppl = init_ppl;  % initial condition must be reset in each run

    for t = 1:t_max
        new_ppl = zeros(1,N);
        for i = 1:N
            [lang_of_point, frac_lang1, frac_lang2] = find_frac(G,ppl,i);
            prob = switching_prob(frac_lang1, frac_lang2,lang_of_point, s, a);
            switches = prob > rand;
            new_ppl(i) = (1-switches)*(lang_of_point) + switches*(1-lang_of_point);
        end
        l2(t,sim) = sum(new_ppl)/N;
        l1(t,sim) = 1 - l2(t,sim);
        ppl = new_ppl;
    end
end
toc

%% GRAPH


for g = 1:sim_max
    plot(1:t_max,l1(:,g))
    hold on;
end
    xlabel('time steps'), ylabel('fraction of lang1');
    legend
    hold off;


% do over many erdos reyni graphs and get vector of fractions over a certain amt
% of time steps
% plot 1000 erdos reyni ones, or plot average/std
% save vector for each erdos reyni
% a few status parameters
% same amt of ppl in both
% make lattice random 0s and 1s

%% FUNCTIONS

% just look at rows, vertices are not described by (i,j) in G!
% update language/people vector, person 1 = row 1 of G
% don't need to save languages, just work on the vector of neighbors
% person 1 neighbors: row 1 of G
% column numbers are different people, if entry is 1 at a column theres a
% connection and need to check thir language in the vector at same value as
% in column
% lang_count (ppl speaking other language) & neighbors (total amt of 1s in that row)

function [G,n,m] = create_ER_Graph(n,p,seed,format)
%% 
%    Description:
%        this function create Erdos-Renyi random Graph*
%    Author: 
%        X.C.
%    Date: 
%        Oct 25 2016
%    Comment:
%        *This code only generate approximately Erdos-Renyi Random Graph. 
%        Since Erdos-Renyi Model only consider the undirected, non-self-loop
%        graphs. However, this code would firstly create a directed graph with,
%        self-loops. And then transform the directed graph into undirected simply
%        by ignore the upper triangular adjacency matrix and delete the self-loops
%        
%        However, when the graph size n is large enough, the generated graph would
%        approximately similar to the expected Erdos-Renyi Model.
%    Output Arguments:
%        G : generated random graph
%        n : graph size, number of vertexes, |V|
%        m : graph size, number of edges, |E|
%    Input Arguments:
%        n : graph size, number of vertexes, |V|
%        p : the probability p of the second definition of Erdos-Renyi model.
%        seed: seed of the function. 
%        format:
%        opt:
% 
switch nargin
    case 2
        seed=0;
        format=1;
        verbose=false;
    case 3
        format=1;
        verbose=false;
    case 4
        verbose=false;
    otherwise
        disp('input argument invalid')
    
end
rng(seed);
G = spones(triu(sprand(n,n,p),1));
if nargout>2
    m = nnz(G);
end
if format==1
    G = G + G';
end
end

function [lang_of_point, frac_lang1, frac_lang2] = find_frac(graph,lang_vector,r)
lang_of_point = lang_vector(r);
valid_neighbors = graph(r,:);  % only works if full matrix is used, otherwise ignores neighbors
% with index less than r (since these are left as 0 in the matrix, in lower
% triangular part)
amt_neighbors = sum(valid_neighbors);

count = sum(valid_neighbors.*lang_vector);

% checks if you have neighbors. If not, sets frac to 0 for both to avoid
% dividing by zero (these points will never switch languages or affect others)
if amt_neighbors > 0
    frac_lang2 = count/amt_neighbors;
    frac_lang1 = 1 - frac_lang2;
else
    frac_lang2 = 0;
    frac_lang1 = 0;
end
end

function [prob_at_point] = switching_prob(frac_lang1, frac_lang2, lang_of_point, s, a)
%P(2 to 1) = sx^a, x = frac_lang1
prob_at_point = (lang_of_point)*s*(frac_lang1)^a + (1-lang_of_point)*(1-s)*(frac_lang2)^a;
end
