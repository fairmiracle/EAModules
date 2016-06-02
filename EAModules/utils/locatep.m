% locate the position in a roulette wheel for selection
function idx = locatep(accumscore,p)
for i = 1:length(accumscore)
    if p < accumscore(i)
        idx = i;
        break;
    end
end