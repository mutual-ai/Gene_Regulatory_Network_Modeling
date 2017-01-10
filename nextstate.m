function nstate = nextstate(con, states,phi)
nstate = (con*states)>phi;
end