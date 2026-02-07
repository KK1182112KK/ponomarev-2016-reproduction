function u_val = interp_control_history(t_query, t_hist, u_hist)
%INTERP_CONTROL_HISTORY  Zero-order hold lookup of u(t_query) from discrete history.
%
%  u_val = interp_control_history(t_query, t_hist, u_hist)
%
%  Inputs:
%    t_query  — scalar time at which u is needed
%    t_hist   — (N x 1) vector of stored time stamps (ascending)
%    u_hist   — (N x 1) vector of stored control values
%
%  Output:
%    u_val    — control value at t_query (zero-order hold: most recent value
%               with time stamp <= t_query)
%
%  If t_query is before the earliest stored time, returns the earliest value.
%  If t_query is after the latest stored time, returns the latest value.

    if isempty(t_hist)
        u_val = 0;
        return;
    end

    if t_query <= t_hist(1)
        u_val = u_hist(1);
        return;
    end

    if t_query >= t_hist(end)
        u_val = u_hist(end);
        return;
    end

    % Binary search for largest index with t_hist(idx) <= t_query
    lo = 1;
    hi = length(t_hist);
    while lo < hi
        mid = lo + floor((hi - lo + 1) / 2);
        if t_hist(mid) <= t_query
            lo = mid;
        else
            hi = mid - 1;
        end
    end
    u_val = u_hist(lo);
end
