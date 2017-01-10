function ostate = overexpress(con, states, genes)
oe = zeros(size(states));
oe(genes,:) = 1; 
states2 = (states + oe)>1e-10; 
ostate = nextstate(con, states2);
end