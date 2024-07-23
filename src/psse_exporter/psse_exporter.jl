_permissive_parse_int(x) = Int64(parse(Float64, x))  # Parses "1.0" as 1, errors on "1.5"

function _increment_gen_char(gen_char::Char)
    (gen_char == '9') && return 'A'
    (gen_char == 'Z') && return '0'
    return gen_char + 1
end

function _increment_gen_id(gen_id::String)
    carry = (last(gen_id) == 'Z')
    if length(gen_id) == 1
        carry && return '0' * _increment_gen_char(first(gen_id))
        return string(_increment_gen_char(first(gen_id)))
    end
    return (carry ? _increment_gen_char(first(gen_id)) : first(gen_id)) *
           _increment_gen_char(last(gen_id))
end

"""
Try to make an informative one or two character name for the generator

  - "generator-1234-AB" -> "AB"
  - "123_CT_7" -> "7"
"""
function _first_choice_gen_id(name::String)
    my_split = argmax(length, [split(name, "-"), split(name, "_")])
    return uppercase(first(last(my_split), 2))
end

"""
Given a vector of Generators, create a corresponding vector of unique-per-bus PSS/E-compatible generator IDs

WRITTEN TO SPEC: PSS/E 33.3 POM page 5-14, PSS/E 34.1 DataFormat page 1-15
"""
function create_gen_ids(gens::Vector{<:PSY.Device})
    gen_ids = Vector{String}()
    sizehint!(gen_ids, length(gens))
    ids_by_bus = Dict{Int64, Vector{String}}()

    for gen in gens
        bus_n = PSY.get_number(PSY.get_bus(gen))
        haskey(ids_by_bus, bus_n) || (ids_by_bus[bus_n] = Vector{String}())
        my_blocked = ids_by_bus[bus_n]
        my_name = _first_choice_gen_id(PSY.get_name(gen))
        while my_name in my_blocked
            my_name = _increment_gen_id(my_name)
        end
        push!(gen_ids, my_name)
        push!(my_blocked, my_name)
    end
    return gen_ids
end

# TODO maybe we want a special System constructor that takes the JSON and does this internally
"""
Given a metadata dictionary parsed from a `raw_metadata_log.json`, yields a function that
can be passed as the `gen_name_formatter` kwarg to the `System` constructor to properly
restore Sienna generator names:
`System("filename.raw"; gen_name_formatter = PF.make_gen_name_formatter_from_metadata(md))`
"""
function make_gen_name_formatter_from_metadata(md::Dict)
    gen_map = md["Gen_Name_Mapping"]
    function md_gen_name_formatter(device_dict::Dict)::String
        bus_n, gen_id = device_dict["source_id"][2:3]
        return gen_map["$(bus_n)_$(rstrip(gen_id))"]
    end
    return md_gen_name_formatter
end

# TODO document this function
function Write_Sienna2PSSE(sys::System, scenario_name::String, year::Int64;
    export_location::Union{Nothing, String} = nothing, base_case = false, setpoint = false,
    setpoint_ts::Union{Nothing, Dates.DateTime} = nothing,
    results_dir::Union{Nothing, String} = nothing, v33 = true)

    # Metadata file
    raw_file_metadata = OrderedDict()

    if (export_location === nothing)
        @warn "Location to save the incremental raw file not specified. Using the data folder in the test folder of the module."
        export_location = joinpath(dirname(dirname(@__DIR__)), "test", "data", "Raw_Export")
    else
        export_location = joinpath(export_location, "Raw_Export")
    end
    @info "Exporting to $export_location"

    if (setpoint)
        if (setpoint_ts === nothing || results_dir === nothing)
            error(
                "The export timestamp and/or sim store results directory must be specifed.",
            )
        end
    end

    raw_exp_dir = joinpath(export_location, scenario_name, string(year))

    # Metadata
    # push!(raw_file_metadata, "Sienna_System_Location" => sys_location)  # TODO fix
    push!(raw_file_metadata, "Raw_File_Export_Location" => raw_exp_dir)
    if base_case
        push!(raw_file_metadata, "Base_Case_Export" => true)
    else
        push!(raw_file_metadata, "Base_Case_Export" => false)
    end

    if ~(isdir(raw_exp_dir))
        mkpath(raw_exp_dir)
    end

    # System should be in NATURAL_UNITS for export (MW)
    PSY.set_units_base_system!(sys, PSY.IS.UnitSystem.NATURAL_UNITS)

    export_ts =
        if (setpoint)
            populate_setpoints!(sys, results_dir, setpoint_ts)
        end

    raw_file = []

    if ~(v33)
        push!(raw_file, "@!IC, SBASE, REV, XFRRAT, NXFRAT, BASFRQ")
    end

    dt_now = Dates.format(Dates.now(), "dd-u-yy-H-M")
    if (v33)
        push!(
            raw_file,
            "0,   $(PSY.get_base_power(sys)), 33, 0, 1, 60.00     /  PSS(R)E 33 RAW created by Sienna2PF.jl  $(dt_now)",
        )
    else
        push!(
            raw_file,
            "0,   $(PSY.get_base_power(sys)), 34, 0, 1, 60.00 /  PSS(R)E 34 RAW created by Sienna2PF.jl  $(dt_now)",
        )
    end

    sys_name =
        if (setpoint)
            "$(scenario_name)_$(string(year))_$(export_ts)"
        else
            "$(scenario_name)_$(string(year))"
        end

    push!(raw_file, "File Name: $(sys_name)")
    push!(raw_file, "         ")

    if ~(v33)
        push!(raw_file, "GENERAL, THRSHZ=0.0001, PQBRAK=0.7, BLOWUP=5.0")
        push!(raw_file, "GAUSS, ITMX=100, ACCP=1.6, ACCQ=1.6, ACCM=1.0, TOL=0.0001")
        push!(
            raw_file,
            "NEWTON, ITMXN=20, ACCN=1.0, TOLN=0.1, VCTOLQ=0.1, VCTOLV=0.00001, DVLIM=0.99, NDVFCT=0.99",
        )
        push!(
            raw_file,
            "ADJUST, ADJTHR=0.005, ACCTAP=1.0, TAPLIM=0.05, SWVBND=100.0, MXTPSS=99, MXSWIM=10",
        )
        push!(raw_file, "TYSL, ITMXTY=20, ACCTY=1.0, TOLTY=0.00001")
        push!(raw_file, """RATING, 1, "RATE1 ", "RATING SET 1                    " """)
        push!(raw_file, """RATING, 2, "RATE2 ", "RATING SET 2                    " """)
        push!(raw_file, """RATING, 3, "RATE3 ", "RATING SET 3                    " """)
        push!(raw_file, "0 / END OF SYSTEM-WIDE DATA, BEGIN BUS DATA")
        push!(
            raw_file,
            "@!   I,'NAME        ', BASKV, IDE,AREA,ZONE,OWNER, VM,        VA,    NVHI,   NVLO,   EVHI,   EVLO",
        )
    end

    # Buses
    # Find buses in PSY System
    psy_buses = collect(PSY.get_components(PSY.Bus, sys))

    bus_mapping = OrderedDict()
    # Metadata
    push!(raw_file_metadata, "Bus_Name_Mapping" => OrderedDict())

    for bus in psy_buses
        _bus_num = PSY.get_number(bus)
        bus_num = _bus_num
        if (bus_num > 1000000)
            bus_num = bus_num - 900000
            if ~(has_bus(sys, bus_num))
                push!(bus_mapping, _bus_num => bus_num)
            else
                while (has_bus(sys, bus_num))
                    bus_num = rand(900000:1000000)
                end
                push!(bus_mapping, _bus_num => bus_num)
            end
            @warn "$(PSY.get_name(bus)) has more than a 6-digit bus number $(_bus_num), changing this to $(bus_num))"
        end

        bus_name = "BUS_$(string(bus_num))"
        push!(raw_file_metadata["Bus_Name_Mapping"], PSY.get_name(bus) => bus_name)
        bus_base_kv = PSY.get_base_voltage(bus)
        bus_ide = PSY.get_bustype(bus).value
        area = _permissive_parse_int(PSY.get_name(PSY.get_area(bus)))
        l_z = _permissive_parse_int(PSY.get_name(PSY.get_load_zone(bus)))
        owner = 1 # DEFAULT
        v_mag = PSY.get_magnitude(bus)
        v_ang = rad2deg(PSY.get_angle(bus))
        #v_ang = PSY.get_angle(bus)

        bus_entry = "$(bus_num), '$(bus_name)', $(bus_base_kv), $(bus_ide), $(area), $(l_z), $(owner), $(v_mag), $(v_ang), 1.1, 0.9, 1.1, 0.9"
        push!(raw_file, bus_entry)
    end

    # Metadata
    push!(raw_file_metadata, "Bus_Number_Mapping" => bus_mapping)

    if (v33)
        push!(raw_file, "0 /End of Bus data, Begin Load data")
    else
        push!(raw_file, "0 / END OF BUS DATA, BEGIN LOAD DATA")
        push!(
            raw_file,
            "@!   I,'ID',STAT,AREA,ZONE,      PL,        QL,        IP,        IQ,        YP,        YQ, OWNER,SCALE,INTRPT,  DGENP,     DGENQ, DGENF",
        )
    end

    # Loads
    # Find buses in PSY System where LOADS are available
    psy_loads = collect(PSY.get_components(PSY.StandardLoad, sys))

    # V34
    dg_enp = 0.0 # DEFAULT
    dg_enq = 0.0 # DEFAULT
    dg_enm = 0 # DEFAULT

    for load in psy_loads
        load_bus_num = PSY.get_number(PSY.get_bus(load))
        if (haskey(bus_mapping, load_bus_num))
            load_bus_num = bus_mapping[load_bus_num]
        end
        load_id = last(PSY.get_name(load))
        load_status = get_PSSE_status(PSY.get_available(load))
        load_bus = PSY.get_bus(load)
        area = _permissive_parse_int(PSY.get_name(PSY.get_area(load_bus)))
        l_z = _permissive_parse_int(PSY.get_name(PSY.get_load_zone(load_bus)))
        p_l = base_case ? 0.0 : PSY.get_constant_active_power(load)
        q_l = base_case ? 0.0 : PSY.get_constant_reactive_power(load)
        PSY.set_units_base_system!(sys, PSY.IS.UnitSystem.DEVICE_BASE) # These need to be on device base to be parsed correctly for now, everything else, sys base
        i_p = base_case ? 0.0 : PSY.get_current_active_power(load)
        i_q = base_case ? 0.0 : PSY.get_current_reactive_power(load)
        y_p = base_case ? 0.0 : PSY.get_impedance_active_power(load)
        y_q = base_case ? 0.0 : PSY.get_impedance_reactive_power(load)
        PSY.set_units_base_system!(sys, PSY.IS.UnitSystem.NATURAL_UNITS) # Back to sys base
        owner = 1 # DEFAULT
        scale = 1 # DEFAULT
        intprt = 0 # DEFAULT

        load_entry =
            if (v33)
                "$(load_bus_num), '$(load_id) ', $(load_status), $(area), $(l_z), $(p_l), $(q_l), $(i_p), $(i_q), $(y_p), $(y_q), $(owner), $(scale), $(intprt)"
            else
                "$(load_bus_num), '$(load_id) ', $(load_status), $(area), $(l_z), $(p_l), $(q_l), $(i_p), $(i_q), $(y_p), $(y_q), $(owner), $(scale), $(intprt), $(dg_enp), $(dg_enq), $(dg_enm)"
            end

        push!(raw_file, load_entry)
    end

    if (v33)
        push!(raw_file, "0 /End of Load data, Begin Fixed shunt data")
    else
        push!(raw_file, "0 / END OF LOAD DATA, BEGIN FIXED SHUNT DATA")
        push!(raw_file, "@!   I, 'ID', STATUS, GL, BL ")
    end

    # Shunts
    psy_shunts = collect(PSY.get_components(PSY.FixedAdmittance, sys))

    for shunt in psy_shunts
        i = PSY.get_number(PSY.get_bus(shunt))
        id = uppercase(first(last(split(PSY.get_name(shunt), "-")), 2))
        stat = get_PSSE_status(PSY.get_available(shunt))
        gl = real(PSY.get_Y(shunt)) * PSY.get_base_power(sys) # Sienna expects system base, but PSS/e expects MW
        bl = imag(PSY.get_Y(shunt)) * PSY.get_base_power(sys) # Sienna expects system base, but PSS/e expects MW

        if (v33)
            shunt_entry = "$(i), $(id), $(stat), $(gl), $(bl)"
        else # FIX THIS FOR 34
            shunt_entry = "$(i), $(id), $(stat), $(gl), $(bl)"
        end

        push!(raw_file, shunt_entry)
    end

    if (v33)
        push!(raw_file, "0 /End of Fixed shunt data, Begin Generator data")
    else
        push!(raw_file, "0 / END OF FIXED SHUNT DATA, BEGIN GENERATOR DATA")
        push!(
            raw_file,
            "@!   I,'ID',      PG,        QG,        QT,        QB,     VS,    IREG,     MBASE,     ZR,         ZX,         RT,         XT,     GTAP,STAT, RMPCT,      PT,        PB,    O1,    F1,  O2,    F2,  O3,    F3,  O4,    F4,WMOD,  WPF",
        )
    end
    # Generators
    psy_gens = collect(PSY.get_components(PSY.Generator, sys))
    sources = collect(PSY.get_components(PSY.Source, sys))

    gen_ids = create_gen_ids(vcat(psy_gens, sources))
    gen_mapping = OrderedDict{String, String}()  # Maps "$(psse_bus_number)_$(psse_gen_id)" to original PSY name

    for (gen, gen_id) in zip(psy_gens, first(gen_ids, length(psy_gens)))
        if gen isa PSY.ThermalStandard
            gen_bus_num = PSY.get_number(PSY.get_bus(gen))
            if (haskey(bus_mapping, gen_bus_num))
                gen_bus_num = bus_mapping[gen_bus_num]
            end
            p_g = base_case ? 0.0 : PSY.get_active_power(gen)
            q_g = base_case ? 0.0 : PSY.get_reactive_power(gen)
            if PSY.get_reactive_power_limits(gen)[:max] > PSY.get_rating(gen)
                q_t = 999.0
                q_b = -999.0
            else
                q_t = PSY.get_reactive_power_limits(gen)[:max]
                q_b = PSY.get_reactive_power_limits(gen)[:min]
            end
            v_s = PSY.get_magnitude(PSY.get_bus(gen))
            ireg = 0 # DEFAULT 
            mbase = PSY.get_base_power(gen)
            z_r = 0.0 # DEFAULT
            z_x = 1.0 # DEFAULT
            r_t = 0.0 # DEFAULT
            x_t = 0.0 # DEFAULT
            gtap = 1.0 # DEFAULT
            stat = get_PSSE_status(PSY.get_available(gen))
            rmpct = 100.0 # DEFAULT
            p_t = PSY.get_active_power_limits(gen)[:max]
            p_b = PSY.get_active_power_limits(gen)[:min]
            o_i = 1 # DEFAULT
            f_i = 1.0 # DEFAULT
            wmod = get_psse_gen_wmod(gen)
            wpf = 1.0 # DEFAULT
            # v34
            n_reg = 0 # DEFAULT

            gen_entry =
                if (v33)
                    "$(gen_bus_num), '$(gen_id) ', $(p_g), $(q_g), $(q_t), $(q_b), $(v_s), $(ireg), $(mbase), $(z_r), $(z_x), $(r_t), $(x_t), $(gtap), $(stat), $(rmpct), $(p_t), $(p_b), $(o_i), $(f_i), $(wmod), $(wpf)"
                else
                    "$(gen_bus_num), '$(gen_id) ', $(p_g), $(q_g), $(q_t), $(q_b), $(v_s), $(ireg), $(mbase), $(z_r), $(z_x), $(r_t), $(x_t), $(gtap), $(stat), $(rmpct), $(p_t), $(p_b), $(o_i), $(f_i), $(wmod), $(wpf), $(n_reg)"
                end

            push!(raw_file, gen_entry)

        else
            gen_bus_num = PSY.get_number(PSY.get_bus(gen))
            if (haskey(bus_mapping, gen_bus_num))
                gen_bus_num = bus_mapping[gen_bus_num]
            end
            p_g = base_case ? 0.0 : PSY.get_active_power(gen)
            q_g = base_case ? 0.0 : PSY.get_reactive_power(gen)
            q_t = PSY.get_rating(gen)cos(π / 4)
            q_b = 0.0
            v_s = PSY.get_magnitude(PSY.get_bus(gen))
            ireg = 0 # DEFAULT 
            mbase = PSY.get_base_power(gen)
            z_r = 0.0 # DEFAULT
            z_x = 1.0 # DEFAULT
            r_t = 0.0 # DEFAULT
            x_t = 0.0 # DEFAULT
            gtap = 1.0 # DEFAULT
            stat = get_PSSE_status(PSY.get_available(gen))
            rmpct = 100.0 # DEFAULT
            p_t = PSY.get_rating(gen)cos(π / 4)
            p_b = 0.0
            o_i = 1 # DEFAULT
            f_i = 1.0 # DEFAULT
            wmod = get_psse_gen_wmod(gen)
            wpf = 1.0 # DEFAULT
            # v34
            n_reg = 0 # DEFAULT

            gen_entry =
                if (v33)
                    "$(gen_bus_num), '$(gen_id) ', $(p_g), $(q_g), $(q_t), $(q_b), $(v_s), $(ireg), $(mbase), $(z_r), $(z_x), $(r_t), $(x_t), $(gtap), $(stat), $(rmpct), $(p_t), $(p_b), $(o_i), $(f_i), $(wmod), $(wpf)"
                else
                    "$(gen_bus_num), '$(gen_id) ', $(p_g), $(q_g), $(q_t), $(q_b), $(v_s), $(ireg), $(mbase), $(z_r), $(z_x), $(r_t), $(x_t), $(gtap), $(stat), $(rmpct), $(p_t), $(p_b), $(o_i), $(f_i), $(wmod), $(wpf), $(n_reg)"
                end
            push!(raw_file, gen_entry)
        end

        gen_mapping["$(gen_bus_num)_$(gen_id)"] = PSY.get_name(gen)
    end
    push!(raw_file_metadata, "Gen_Name_Mapping" => gen_mapping)

    for (source, gen_id) in zip(sources, last(gen_ids, length(sources)))
        gen_bus_num = PSY.get_number(PSY.get_bus(source))
        if (haskey(bus_mapping, gen_bus_num))
            gen_bus_num = bus_mapping[gen_bus_num]
        end
        p_g = base_case ? 0.0 : PSY.get_active_power(source) * PSY.get_base_power(sys)
        q_g = base_case ? 0.0 : PSY.get_reactive_power(source) * PSY.get_base_power(sys)
        q_t = 10000 * cos(π / 4)
        q_b = 0.0
        v_s = PSY.get_internal_voltage(source)
        ireg = 0 # DEFAULT 
        mbase = PSY.get_base_power(source)
        z_r = 0.0 # DEFAULT
        z_x = 1.0 # DEFAULT
        r_t = 0.0 # DEFAULT
        x_t = 0.0 # DEFAULT
        gtap = 1.0 # DEFAULT
        stat = get_PSSE_status(PSY.get_available(source))
        rmpct = 100.0 # DEFAULT
        p_t = 10000 * cos(π / 4)
        p_b = 0.0
        o_i = 1 # DEFAULT
        f_i = 1.0 # DEFAULT
        wmod = 0
        wpf = 1.0 # DEFAULT
        # v34
        n_reg = 0 # DEFAULT

        gen_entry =
            if (v33)
                "$(gen_bus_num), '$(gen_id) ', $(p_g), $(q_g), $(q_t), $(q_b), $(v_s), $(ireg), $(mbase), $(z_r), $(z_x), $(r_t), $(x_t), $(gtap), $(stat), $(rmpct), $(p_t), $(p_b), $(o_i), $(f_i), $(wmod), $(wpf)"
            else
                "$(gen_bus_num), '$(gen_id) ', $(p_g), $(q_g), $(q_t), $(q_b), $(v_s), $(ireg), $(mbase), $(z_r), $(z_x), $(r_t), $(x_t), $(gtap), $(stat), $(rmpct), $(p_t), $(p_b), $(o_i), $(f_i), $(wmod), $(wpf), $(n_reg)"
            end
        push!(raw_file, gen_entry)

        gen_mapping["$(gen_bus_num)_$(gen_id)"] = PSY.get_name(source)
    end

    if (v33)
        push!(raw_file, "0 /End of Generator data, Begin Branch data")
    else
        push!(raw_file, "0 / END OF GENERATOR DATA, BEGIN BRANCH DATA")
        push!(
            raw_file,
            "@!   I,     J,'CKT',     R,          X,         B,                    'N A M E'                 ,   RATE1,   RATE2,   RATE3,    GI,       BI,       GJ,       BJ,STAT,MET,  LEN,  O1,  F1,    O2,  F2,    O3,  F3,    O4,  F4",
        )
    end

    # Branch
    psy_branches = collect(
        PSY.get_components(
            (x -> PSY.get_from_bus(x) ∈ psy_buses || PSY.get_to_bus(x) ∈ psy_buses),
            PSY.ACBranch,
            sys,
        ),
    )

    # line_dict = Dict{Tuple{Int64, Int64}, Int64}()
    for branch in psy_branches
        if branch isa PSY.Line
            i = PSY.get_number(PSY.get_from_bus(branch))
            if (haskey(bus_mapping, i))
                i = bus_mapping[i]
            end

            j = PSY.get_number(PSY.get_to_bus(branch))
            if (haskey(bus_mapping, j))
                j = bus_mapping[j]
            end
            # if haskey(line_dict, (i,j)) == false
            #     line_dict[(i, j)] = 1
            # elseif haskey(line_dict, (i,j)) == true
            #     line_dict[(i,j)] = (get(line_dict, (i,j), missing)+ 1)
            # end
            # ckt = get(line_dict, (i,j), missing)
            ckt = last(split(PSY.get_name(branch), "_"))
            r = PSY.get_r(branch)
            x = PSY.get_x(branch)
            b = getfield(PSY.get_b(branch), :from) + getfield(PSY.get_b(branch), :to)
            rate_a = PSY.get_rating(branch)
            rate_b = PSY.get_rating(branch)
            rate_c = PSY.get_rating(branch)
            g_i, b_i = 0.0, 0.0 # DEFAULT
            g_j, b_j = 0.0, 0.0 # DEFAULT
            st = get_PSSE_status(PSY.get_available(branch))
            met = 1 # DEFAULT
            len = 0.0  # DEFAULT
            o_i = 1 # DEFAULT
            f_i = 1.0 # DEFAULT

            # V34
            name = filter(
                x -> !isspace(x),
                replace(first(PSY.get_name(branch), 40), r"," => ""),
            )

            branch_entry =
                if (v33)
                    "$(i), $(j), '$(ckt)', $(r), $(x), $(b), $(rate_a), $(rate_b), $(rate_c), $(g_i), $(b_i), $(g_j), $(b_j), $(st), $(met), $(len), $(o_i), $(f_i)"
                else
                    "$(i), $(j), '$(ckt)', $(r), $(x), $(b), $(name), $(rate_a), $(rate_b), $(rate_c), $(g_i), $(b_i), $(g_j), $(b_j), $(st), $(met), $(len), $(o_i), $(f_i)"
                end

            push!(raw_file, branch_entry)
        else
            continue
        end
    end

    # Metadata
    push!(raw_file_metadata, "Other_ID_Legend_Mapping" => Dict())
    push!(raw_file_metadata["Other_ID_Legend_Mapping"], "ENG" => "ENGAGE_Bus")
    push!(raw_file_metadata["Other_ID_Legend_Mapping"], "EL" => "ENGAGE_Line")

    if (v33)
        push!(raw_file, "0 /End of Branch data, Begin Transformer data")
    else
        push!(raw_file, "0 / END OF BRANCH DATA, BEGIN SYSTEM SWITCHING DEVICE DATA")
        push!(raw_file, "0 / END OF SYSTEM SWITCHING DEVICE DATA, BEGIN TRANSFORMER DATA")
        push!(
            raw_file,
            "@! I, J, K, CKT, CW, CZ, CM, MAG1, MAG2, NMETER, 'NAME', STAT, O1, F1, O2, F2, O3, F3, O4, F4, VECGRP",
        )
        push!(raw_file, "@! R1-2, X1-2, SBASE1-2")
        push!(
            raw_file,
            "@! WINDV1, NOMV1, ANG1, RATE11, RATE12, RATE13, COD1, CONT1, RMA1, RMI1, VMA1, VMI1",
        )
        push!(raw_file, "@! WINDV2, NOM2")
    end

    # Transformers
    # twowind_xfmr_dict = Dict{Tuple{Int64, Int64}, Int64}()
    # tap_xfmr_dict = Dict{Tuple{Int64, Int64}, Int64}()
    for branch in psy_branches
        if branch isa PSY.Transformer2W
            i = PSY.get_number(PSY.get_from_bus(branch))
            if (haskey(bus_mapping, i))
                i = bus_mapping[i]
            end

            j = PSY.get_number(PSY.get_to_bus(branch))
            if (haskey(bus_mapping, j))
                j = bus_mapping[j]
            end

            k = 0

            # if haskey(twowind_xfmr_dict, (i,j)) == false
            #     twowind_xfmr_dict[(i, j)] = 1
            # elseif haskey(twowind_xfmr_dict, (i,j)) == true
            #     twowind_xfmr_dict[(i,j)] = (get(twowind_xfmr_dict, (i,j), missing)+ 1)
            # end
            #ckt = get(twowind_xfmr_dict, (i,j), missing)
            ckt = last(split(PSY.get_name(branch), "_"))
            cw = 1 # DEFAULT 
            cz = 1 # DEFAULT 
            cm = 1 # DEFAULT 
            mag1 = 0.0 # DEFAULT  - can calculate?
            mag2 = 0.0 # DEFAULT
            n_meter = 2 # DEFAULT
            name = "$(i)-$(j)_$(ckt)"
            stat = get_PSSE_status(PSY.get_available(branch))
            o1 = 1 # DEFAULT
            f1 = 1 # DEFAULT
            o2 = 0 # DEFAULT
            f2 = 0 # DEFAULT
            o3 = 0 # DEFAULT
            f3 = 0 # DEFAULT
            o4 = 0 # DEFAULT
            f4 = 0 # DEFAULT
            vecgrp = "            " # DEFAULT (12 blanks)
            r1_2 = PSY.get_r(branch)
            x1_2 = PSY.get_x(branch)
            sbase1_2 = PSY.get_base_power(branch)
            windv1 = 1.0 # DEFAULT
            nomv1 = 0.0 # DEFAULT
            ang1 = 0.0 # DEFAULT
            rata1 = PSY.get_rating(branch)
            ratb1 = PSY.get_rating(branch)
            ratc1 = PSY.get_rating(branch)
            cod1 = 0 # DEFAULT
            cont1 = 0 # DEFAULT
            rma1 = 1.1 # DEFAULT
            rmi1 = 0.9 # DEFAULT
            vma1 = 1.1 # DEFAULT
            vmi1 = 0.9 # DEFAULT
            ntp1 = 33 # DEFAULT
            tab1 = 0  # DEFAULT
            cr1 = 0.0 # DEFAULT
            cx1 = 0.0 # DEFAULT
            cnxa1 = 0.0  # DEFAULT
            windvu2 = 1.0 # DEFAULT
            nomv2 = 0.0 # DEFAULT
            #v34
            node1 = 0 # DEFAULT

            branch_entry1 = "$(i), $(j), $(k), $(ckt), $(cw), $(cz), $(cm), $(mag1), $(mag2), $(n_meter), $(name), $(stat), $(o1), $(f1), $(o2), $(f2), $(o3), $(f3) , $(o4), $(f4), $(vecgrp)"

            branch_entry2 = "$(r1_2), $(x1_2), $(sbase1_2)"

            branch_entry3 =
                if (v33)
                    "$(windv1), $(nomv1), $(ang1), $(rata1), $(ratb1), $(ratc1), $(cod1), $(cont1), $(rma1), $(rmi1), $(vma1), $(vmi1), $(ntp1), $(tab1), $(cr1), $(cx1), $(cnxa1)"
                else
                    "$(windv1), $(nomv1), $(ang1), $(rata1), $(ratb1), $(ratc1), $(cod1), $(cont1), $(rma1), $(rmi1), $(vma1), $(vmi1), $(ntp1), $(tab1), $(cr1), $(cx1), $(cnxa1), $(node1)"
                end
            branch_entry4 = "$(windvu2), $(nomv2)"

            push!(raw_file, branch_entry1)
            push!(raw_file, branch_entry2)
            push!(raw_file, branch_entry3)
            push!(raw_file, branch_entry4)

        elseif branch isa PSY.TapTransformer
            i = PSY.get_number(PSY.get_from_bus(branch))
            if (haskey(bus_mapping, i))
                i = bus_mapping[i]
            end

            j = PSY.get_number(PSY.get_to_bus(branch))
            if (haskey(bus_mapping, j))
                j = bus_mapping[j]
            end

            k = 0

            # if haskey(tap_xfmr_dict, (i,j)) == false
            #     tap_xfmr_dict[(i, j)] = 1
            # elseif haskey(tap_xfmr_dict, (i,j)) == true
            #     tap_xfmr_dict[(i,j)] = (get(tap_xfmr_dict, (i,j), missing)+ 1)
            # end
            #ckt = get(tap_xfmr_dict, (i,j), missing)
            ckt = last(split(PSY.get_name(branch), "_"))
            cw = 1 # DEFAULT 
            cz = 1 # DEFAULT 
            cm = 1 # DEFAULT 
            mag1 = PSY.get_primary_shunt(branch)#0.0 # DEFAULT  - can calculate?
            mag2 = 0.0 # DEFAULT
            n_meter = 2 # DEFAULT
            name = PSY.get_name(branch)
            stat = get_PSSE_status(PSY.get_available(branch))
            o1 = 1 # DEFAULT
            f1 = 1 # DEFAULT
            o2 = 0 # DEFAULT
            f2 = 0 # DEFAULT
            o3 = 0 # DEFAULT
            f3 = 0 # DEFAULT
            o4 = 0 # DEFAULT
            f4 = 0 # DEFAULT
            vecgrp = "            " # DEFAULT (12 blanks)
            r1_2 = PSY.get_r(branch)
            x1_2 = PSY.get_x(branch)
            sbase1_2 = PSY.get_base_power(branch)
            windv1 = 1.0 # DEFAULT
            nomv1 = 0.0 # DEFAULT
            ang1 = 0.0 # DEFAULT
            rata1 = PSY.get_rating(branch)
            ratb1 = PSY.get_rating(branch)
            ratc1 = PSY.get_rating(branch)
            cod1 = 1
            cont1 = j
            rma1 = 1.1 # DEFAULT
            rmi1 = 0.9 # DEFAULT
            vma1 = 1.1 # DEFAULT
            vmi1 = 0.9 # DEFAULT
            ntp1 = 33 # DEFAULT
            tab1 = 0  # DEFAULT
            cr1 = 0.0 # DEFAULT
            cx1 = 0.0 # DEFAULT
            cnxa1 = 0.0  # DEFAULT
            windvu2 = 0.99
            nomv2 = 0.0 # DEFAULT

            branch_entry1 = "$(i), $(j), $(k), $(ckt), $(cw), $(cz), $(cm), $(mag1), $(mag2), $(n_meter), $(name), $(stat), $(o1), $(f1), $(o2), $(f2), $(o3), $(f3) , $(o4), $(f4), $(vecgrp)"
            branch_entry2 = "$(r1_2), $(x1_2), $(sbase1_2)"
            branch_entry3 =
                if (v33)
                    "$(windv1), $(nomv1), $(ang1), $(rata1), $(ratb1), $(ratc1), $(cod1), $(cont1), $(rma1), $(rmi1), $(vma1), $(vmi1), $(ntp1), $(tab1), $(cr1), $(cx1), $(cnxa1)"
                else
                    "$(windv1), $(nomv1), $(ang1), $(rata1), $(ratb1), $(ratc1), $(cod1), $(cont1), $(rma1), $(rmi1), $(vma1), $(vmi1), $(ntp1), $(tab1), $(cr1), $(cx1), $(cnxa1), $(node1)"
                end
            branch_entry4 = "$(windvu2), $(nomv2)"

            push!(raw_file, branch_entry1)
            push!(raw_file, branch_entry2)
            push!(raw_file, branch_entry3)
            push!(raw_file, branch_entry4)
        end
    end

    if (v33)
        push!(raw_file, "0 /End of Transformer data, Begin Area interchange data")
        push!(raw_file, "0 /End of Area interchange data, Begin Two-terminal dc line data")
        push!(raw_file, "0 /End of Two-terminal dc line data, Begin VSC dc line data")
        push!(raw_file, "0 /End of VSC dc line data, Begin Impedance correction table data")
        push!(
            raw_file,
            "0 /End of Impedance correction table data, Begin Multi-terminal dc line data",
        )
        push!(
            raw_file,
            "0 /End of Multi-terminal dc line data, Begin Multi-section line data",
        )
        push!(raw_file, "0 /End of Multi-section line data, Begin Zone data")
        push!(raw_file, "0 /End of Zone data, Begin Inter-area transfer data")
        push!(raw_file, "0 /End of Inter-area transfer data, Begin Owner data")
        push!(raw_file, "0 /End of Owner data, Begin FACTS device data")
        push!(raw_file, "0 /End of FACTS device data, Begin Switched shunt data")
        push!(raw_file, "0 /End of Switched shunt data, Begin GNE device data")
        push!(raw_file, "0 /End of GNE device data, Begin Induction machine data")
        push!(raw_file, "0 /End of Induction machine data")
    else
        push!(raw_file, "0 / END OF TRANSFORMER DATA, BEGIN AREA DATA")
        push!(raw_file, "0 / END OF AREA DATA, BEGIN TWO-TERMINAL DC DATA")
        push!(raw_file, "0 / END OF TWO-TERMINAL DC DATA, BEGIN VSC DC LINE DATA")
        push!(raw_file, "0 / END OF VSC DC LINE DATA, BEGIN IMPEDANCE CORRECTION DATA")
        push!(
            raw_file,
            "0 / END OF IMPEDANCE CORRECTION DATA, BEGIN MULTI-TERMINAL DC DATA",
        )
        push!(raw_file, "0 / END OF MULTI-TERMINAL DC DATA, BEGIN MULTI-SECTION LINE DATA")
        push!(raw_file, "0 / END OF MULTI-SECTION LINE DATA, BEGIN ZONE DATA")
        push!(raw_file, "0 / END OF ZONE DATA, BEGIN INTER-AREA TRANSFER DATA")
        push!(raw_file, "0 / END OF INTER-AREA TRANSFER DATA, BEGIN OWNER DATA")
        push!(raw_file, "0 / END OF OWNER DATA, BEGIN FACTS DEVICE DATA")
        push!(raw_file, "0 / END OF FACTS DEVICE DATA, BEGIN SWITCHED SHUNT DATA")
        push!(raw_file, "0 / END OF SWITCHED SHUNT DATA, BEGIN GNE DATA")
        push!(raw_file, "0 / END OF GNE DATA, BEGIN INDUCTION MACHINE DATA")
        push!(raw_file, "0 / END OF INDUCTION MACHINE DATA, BEGIN SUBSTATION DATA")
        push!(raw_file, "0 / END OF SUBSTATION DATA")
    end

    push!(raw_file, "Q")

    @info "Exporting raw file and relevant metadata log ..."

    # Export raw file
    export_raw_file_loc =
        if (setpoint)
            joinpath(raw_exp_dir, "$(scenario_name)_$(export_ts).raw")
        else
            joinpath(raw_exp_dir, "$(scenario_name).raw")
        end

    open(export_raw_file_loc, "w") do file
        DelimitedFiles.writedlm(file, raw_file)
    end

    # Exporting the metadata file
    raw_log_json_location =
        if (setpoint)
            joinpath(raw_exp_dir, "raw_metadata_log_$(export_ts).json")
        else
            joinpath(raw_exp_dir, "raw_metadata_log.json")
        end

    open(raw_log_json_location, "w") do f
        JSON.print(f, raw_file_metadata, 4)
    end

    @info "Sienna => PSSE(R) Done."
end
