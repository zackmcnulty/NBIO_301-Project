function r = activation_Acurrent(V)

r = 0.1*(V + 65)/(1 - exp(-(V + 65)/10));


end

