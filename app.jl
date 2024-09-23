module App
# == Packages ==
# set up Genie development environment. Use the Package Manager to install new packages
using GenieFramework, DataFrames, PlotlyBase
@genietools

restriction_enzymes = Dict(
    "EcoRI" => (site = "GAATTC", cut_pos = 0),
    "HindIII" => (site = "AAGCTT", cut_pos = 0),
    "BamHI" => (site = "GGATCC", cut_pos = 0),
    "PstI" => (site = "CTGCAG", cut_pos = 4),
    "NotI" => (site = "GCGGCCGC", cut_pos = 1),
    "SalI" => (site = "GTCGAC", cut_pos = 0),
    "EcoRV" => (site = "GATATC", cut_pos = 0),
    "BglII" => (site = "AGATCT", cut_pos = 0),
    "SmaI" => (site = "CCCGGG", cut_pos = 2),
    "XhoI" => (site = "CTCGAG", cut_pos = 0),
    "KpnI" => (site = "GGTACC", cut_pos = 4),
    "AluI" => (site = "AGCT", cut_pos = 1),
    "XbaI" => (site = "TCTAGA", cut_pos = 0))


function find_cut_sites(dna_seq::String, enzyme::String)
    site, cut_pos = restriction_enzymes[enzyme]
    positions = findall(site, dna_seq)
    cut_sites = [first(pos) + cut_pos for pos in positions]
    return cut_sites, positions
end

# Función para cortar la secuencia de ADN en los sitios de corte
function cut_sequence(dna_seq::String, cut_sites::Vector{Int})
    cut_sites = sort(cut_sites)
    fragments = Vector{String}()
    start_pos = 1
 for cut_site in cut_sites
        push!(fragments, dna_seq[start_pos:cut_site])
        start_pos = cut_site + 1
    end
    push!(fragments, dna_seq[start_pos:end])
    return fragments
end

function generate_marker_bands()
    return [15000, 10000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 850, 650, 500, 400,300, 200, 100]
end


# == Reactive code ==
# add reactive code to make the UI interactive
@app begin
    # == Reactive variables ==
    # reactive variables exist in both the Julia backend and the browser with two-way synchronization
    # @out variables can only be modified by the backend
    # @in variables can be modified by both the backend and the browser
    # variables must be initialized with constant values, or variables defined outside of the @app block
    @in dna_seq = ""
    @in odna_seq = ""
    @in enz1 = ""
    @in enz2 = ""
    @in enz3 = ""
    @out enzymes = ["", "AluI", "BamHI", "BglII", "EcoRI", "EcoRV", "HindIII", "KpnI", "NotI", "PstI", "SalI", "SmaI", "XbaI", "XhoI"]
    @out enzopt2 = String[]
    @out enz3opt = String[]
    @out fragments = String[]
    @out fragments_len = Int[]
    @out total_len::Int64 = 0
    @out table = DataTable(DataFrame(counter=Int[], fragments = String[], fragments_len = Int[]))
    @out marker_bands=[]
    @out sample_bands=[]
    @out x1=[]
    @out x2=[]
    @out marker_colors = ["gray", "gray", "gray", "Gray", "Gray", "Gray", "Gray", "Gray", "Gray", "White", "Gray", "Gray", "Gray", "Gray", "gray", "gray", "gray", "gray"]
    @out marker_line_widths = [2, 2, 2, 2, 3, 3, 3, 3, 4, 5, 4, 3, 3, 4, 5, 5, 6, 6]
    #@out trace1=[]
    #@out trace2=[]
    #@out plotdata = [trace1, trace2]
    #@out plotlayout = layout

    # == Reactive handlers ==
    # reactive handlers watch a variable and execute a block of code when
    # its value changes
    @onchange odna_seq, enz1, enz2, enz3 begin

        dna_seq = replace(uppercase(odna_seq), r"\s+" => "")
        total_len = length(dna_seq)
        selected_enzymes = filter(!isempty, [enz1, enz2, enz3])
        # Validación para asegurar que enz2 es diferente de enz1 y enz3 es diferente de enz1 y enz2
        enzopt2 = filter(e -> e != enz1, enzymes)

        enz3opt = filter(e -> e != enz2, enzopt2)

        #if enz2 == enz1 || (!isempty(enz3) && (enz3 == enz1 || enz3 == enz2))
        #    error("enzymes must be diferent between eachother.")
        #end
        all_cut_sites = Int[]
        enzyme_positions = Dict{String, Vector{Union{Missing, Int}}}()

            # Calcular los fragmentos y posiciones de corte
        for enzyme in selected_enzymes
            cut_sites, positions = find_cut_sites(dna_seq, enzyme)
            append!(all_cut_sites, cut_sites)
            all_cut_sites = unique(sort(all_cut_sites))  # Asegura sitios de corte únicos y ordenados
            fragments = cut_sequence(dna_seq, all_cut_sites)

            # Asegurar que el vector de posiciones tenga la misma longitud que el número de fragmentos
            positions_vector = Vector{Union{Missing, Int}}(missing, length(fragments))
            for (i, pos) in enumerate(positions)
                if i <= length(positions_vector)
                    positions_vector[i] = first(pos)
                end
            end

            enzyme_positions[enzyme] = positions_vector
        end

        # Calcular longitudes de los fragmentos
        fragments_len = [length(frag) for frag in fragments]
        fragments_len = sort(fragments_len)
        fragments = sort(fragments, by = length)

        # Crear el DataFrame y añadir posiciones de cada enzima
        df = DataFrame(No = 1:length(fragments), Fragments = fragments, Length = fragments_len) 

        # Asegurar que las columnas de posiciones tengan la misma longitud que las otras columnas
        for enzyme in selected_enzymes
            pos_vector = enzyme_positions[enzyme]
            if length(pos_vector) < length(fragments)
                append!(pos_vector, fill(missing, length(fragments) - length(pos_vector)))
            end
            # Convertir los valores `missing` a "--"
            df[!, enzyme] = coalesce.(pos_vector, "")
        end

        table = DataTable(df)

        length_vector = df.Length
        # Función para generar las posiciones de las bandas en la muestra digerida
        function generate_sample_bands()
            return length_vector  # Supongamos que estos son los fragmentos resultantes de la digestión
        end

        # Crear DataFrames para las bandas
        marker_df = DataFrame(x = fill(1, length(generate_marker_bands())), y = generate_marker_bands())
        sample_df = DataFrame(x = fill(1.5, length(generate_sample_bands())), y = generate_sample_bands())
        # Función para crear el gráfico del gel de electroforesis
        

        marker_bands = generate_marker_bands()
        sample_bands = generate_sample_bands()
        x1 = marker_df.x
        x2 = sample_df.x
                

    end
    # the onbutton handler will set the variable to false after the block is executed
    
end



# == Pages ==
# register a new route and the page that will be loaded on access
@page("/", "app.jl.html")
end

# == Advanced features ==
#=
- The @private macro defines a reactive variable that is not sent to the browser. 
This is useful for storing data that is unique to each user session but is not needed
in the UI.
    @private table = DataFrame(a = 1:10, b = 10:19, c = 20:29)

=#
