function [f,g,H] = zfun(zk, tk, mc, sc, o)
        a = exp(zk).*o;
        zc = zk - mc;
        sci = 1/sc;
        
        f = tk.*zk - a - 0.5*sci*zc.^2;
        
        if nargout > 1
            g = tk - a - zc*sci;
            
            if nargout > 2
                H = -a - sci;
                %H = spdiags(h, 0, length(h), length(h));
            end
        end
        
    end