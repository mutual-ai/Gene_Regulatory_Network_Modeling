function ostates = overexpress_final(con, states, genes, expresstime)
os_prev = states;
os = overexpress(con,os_prev,genes);
n = 1; 
while ~isequal(os, os_prev) && n<expresstime
    os_prev = os;
    os = overexpress(con,os_prev,genes);
    n = n+1;
end
os_prev = os;
os = nextstate(con,os_prev);
n = 1; maxn = 100;
while ~isequal(os, os_prev) && n<maxn
    os_prev = os;
    os = nextstate(con,os_prev);
    n = n+1;
end
ostates = os;
end