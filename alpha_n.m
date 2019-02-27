function r = alpha_n(V)
    % calculates value of gating function; parameters based on experimental
    % data by Hodgkin-Huxley

    r = 0.01*(V + 55) / (1-exp(-(V+55)/10));
end

