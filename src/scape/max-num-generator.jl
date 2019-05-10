using Distributed
function producer(c::Channel)
    counter = 999999
    while true
        put!(c, counter)
        counter += 1
    end        
end;

