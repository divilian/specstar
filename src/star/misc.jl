
# Misc helper functions.

using Gadfly
using ColorSchemes, Colors

# verbosity: 0=critical, 1=informational, 2=debug
# "prc": "print critical"
# "pri": "print informational"
# "prd": "print debug"
global verbosity = 0

function prc(x) 
    if verbosity ≥ 0
        print(x)
    end
end

function pri(x) 
    if verbosity > 0
        print(x)
    end
end

function prd(x) 
    if verbosity > 1
        print(x)
    end
end

theme = Theme(default_color=colorant"navy")
