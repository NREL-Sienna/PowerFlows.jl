function is_available_source(x, bus::PSY.Bus)
    return PSY.get_available(x) && x.bus == bus
end

function is_available_source(x::PSY.ThermalGen, bus::PSY.Bus)
    return PSY.get_available(x) && x.bus == bus && PSY.get_status(x)
end

function is_available_source(x::PSY.ElectricLoad, bus::PSY.Bus)
    return false
end
