function exportToPly(msh::Mesh, fn::String)
    vts = msh.vertices
    fcs = msh.faces
    nV = size(vts,1)
    nF = size(fcs,1)

    str = open(fn,"w")

    # write the header
    write(str,"ply\n")
    write(str,"format binary_little_endian 1.0\n")
    write(str,"element vertex $nV\n")
    write(str,"property float x\nproperty float y\nproperty float z\n")
    write(str,"element face $nF\n")
    write(str,"property list uchar int vertex_index\n")
    write(str,"end_header\n")

    # write the data
    for i = 1:nV
        v = vts[i]
        write(str,float32(v.e1))
        write(str,float32(v.e2))
        write(str,float32(v.e3))
    end

   for i = 1:nF
        f = fcs[i]
        write(str,uint8(3))
        write(str,int32(f.v1-1))
        write(str,int32(f.v2-1))
        write(str,int32(f.v3-1))
    end
    close(str)
end
export exportToPly

function exportToStl(msh::Mesh, fn::String)
    vts = msh.vertices
    fcs = msh.faces
    nV = size(vts,1)
    nF = size(fcs,1)

    str = open(fn,"w")

    # write the header
    write(str,"solid vcg\n")

    # write the data
    for i = 1:nF
        f = fcs[i]
        n = [0,0,0] # TODO: properly compute normal(f)
        txt = @sprintf "  facet normal %e %e %e\n" n[1] n[2] n[3]
        write(str,txt)
        write(str,"    outer loop\n")
        v = vts[f.v1]
        txt = @sprintf "      vertex  %e %e %e\n" v[1] v[2] v[3]
        write(str,txt)

        v = vts[f.v2]
        txt = @sprintf "      vertex  %e %e %e\n" v[1] v[2] v[3]
        write(str,txt)

        v = vts[f.v3]
        txt = @sprintf "      vertex  %e %e %e\n" v[1] v[2] v[3]
        write(str,txt)

        write(str,"    endloop\n")
        write(str,"  endfacet\n")
    end

    write(str,"endsolid vcg\n")
    close(str)
end
export exportToStl

# | Read a .2dm (SMS Aquaveo) mesh-file and construct a @Mesh@
function import2dm(file::String)
    parseNode(w::Array{String}) = Vertex3(float64(w[3]), float64(w[4]), float64(w[5]))
    parseTriangle(w::Array{String}) = Face(int64(w[3]), int64(w[4]), int64(w[5]))
    # Qudrilateral faces are split up into triangles
    function parseQuad(w::Array{String})
        w[7] = w[3]                     # making a circle
        Face[Face(int64(w[i]), int64(w[i+1]), int64(w[i+2])) for i = [3,5]]
    end 
    con = open(file, "r")
    nd =  Array(Vertex3, 0)
    ele = Array(Face,0)
    for line = readlines(con)
        line = chomp(line)
        w = split(line)
        if w[1] == "ND"
            push!(nd, parseNode(w))
        elseif w[1] == "E3T"
            push!(ele,parseTriangle(w))
        elseif w[1] == "E4Q"
            append!(ele, parseQuad(w))
        else
            continue
        end
    end
    close(con)
    Mesh(nd,ele)
end
export import2dm

# | Write @Mesh@ to an IOStream
function exportTo2dm(f::IO,m::Mesh)
    function renderVertex(i::Int,v::Vertex)
        "ND $i $(v.x) $(v.y) $(v.z)\n"
    end
    function renderFace(i::Int, f::Face)
        "E3T $i $(f.v1) $(f.v2) $(f.v3) 0\n"
    end
    for i = 1:length(m.faces)
        write(f, renderFace(i, m.faces[i]))
    end
    for i = 1:length(m.vertices)
        write(f, renderVertex(i, m.vertices[i]))
    end
    nothing
end

# | Write a @Mesh@ to file in SMS-.2dm-file-format
function exportTo2dm(f::String,m::Mesh)
    con = open(f, "w")
    exportTo2dm(con, m)
    close(con)
    nothing
end
export exportTo2dm

## Msh support
## Based on https://github.com/scharris/WGFEA/blob/master/TMesh.jl
## Copyright (c) 2013 Stephen C Harris
## MIT License

typealias ElTypeNum Uint8
eltypenum(i::Integer) = if i>=1 && i<=typemax(ElTypeNum) convert(ElTypeNum, i) else error("invalid element type: $i") end
eltypenum(s::String) = eltypenum(int(s))

# Gmsh triangle type codes
const ELTYPE_3_NODE_TRIANGLE = eltypenum(2)
# lower order elements, which can be ignored
const ELTYPE_POINT = eltypenum(15)
const ELTYPE_2_NODE_LINE = eltypenum(1)
const ELTYPE_3_NODE_LINE = eltypenum(8)
const ELTYPE_4_NODE_LINE = eltypenum(26)
const ELTYPE_5_NODE_LINE = eltypenum(27)
const ELTYPE_6_NODE_LINE = eltypenum(28)

function is_3_node_triangle_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_3_NODE_TRIANGLE
end

function is_lower_order_el_type(el_type::ElTypeNum)
  el_type == ELTYPE_POINT ||
  el_type == ELTYPE_2_NODE_LINE ||
  el_type == ELTYPE_3_NODE_LINE ||
  el_type == ELTYPE_4_NODE_LINE ||
  el_type == ELTYPE_5_NODE_LINE ||
  el_type == ELTYPE_6_NODE_LINE
end

const TOKN_NODELINE_ELNUM = 1
const TOKN_NODELINE_POINT1 = 2
const TOKN_NODELINE_POINT2 = 3
const TOKN_NODELINE_POINT3 = 4

# Gmsh element line format:
# elm-number elm-type number-of-tags < tag > ... vert-number-list
const TOKN_ELLINE_ELTYPE = 2
const TOKN_ELLINE_NUMTAGS = 3

# Gmsh base elements iterator

immutable GmshElementsIter
  first_base_tri_line::String
  is::IO # stream should be positioned after first base element line
end

# Represents an element provided to the mesh constructor which is to be subdivided to form
# some of the final mesh elements.
immutable BaseTri
    point_nums::Face
    tag_physreg::Int64
    tag_geoment::Int64
    other_tags::Array{Int64,1}
end

const EMPTY_OTHER_TAGS = Array(Int64,0)

function base_tri_from_gmsh_el(toks::Array)
    const num_tags = int(toks[TOKN_ELLINE_NUMTAGS])
    assert(num_tags >= 2)
    const physreg_tag = int64(toks[TOKN_ELLINE_NUMTAGS+1])
    const geoment_tag = int64(toks[TOKN_ELLINE_NUMTAGS+2])
    const other_tags = num_tags >= 3 ? map(t -> int64(t), toks[TOKN_ELLINE_NUMTAGS+3:TOKN_ELLINE_NUMTAGS+num_tags]) :
                                     EMPTY_OTHER_TAGS
    const point_nums = let last_tag_tokn = TOKN_ELLINE_NUMTAGS + num_tags;
        Face(parseint(toks[last_tag_tokn + 1]),
         parseint(toks[last_tag_tokn + 2]),
         parseint(toks[last_tag_tokn + 3]))
    end

    BaseTri(point_nums,
          physreg_tag,
          geoment_tag,
          other_tags)
end

import Base.start, Base.done, Base.next
function start(gmsh_els_iter::GmshElementsIter)
    base_tri_from_gmsh_el(split(gmsh_els_iter.first_base_tri_line, ' '))
end

function done(gmsh_els_iter::GmshElementsIter, next_base_tri)
    next_base_tri == nothing
end

function next(gmsh_els_iter::GmshElementsIter, next_base_tri)
    while true
        const line = readline(gmsh_els_iter.is)
        if beginswith(line, "\$EndElements")
            return (next_base_tri, nothing)
        elseif line == ""
            error("No EndElements marker found.")
        else
            const toks = split(line, ' ')
            const el_type = eltypenum(toks[TOKN_ELLINE_ELTYPE])
            if is_lower_order_el_type(el_type)
                continue
            elseif is_3_node_triangle_el_type(el_type)
                return (next_base_tri, base_tri_from_gmsh_el(toks))
            else
                error("Element type $el_type is not supported.")
            end
        end
    end
end

function importmsh(io::IO)
    base_pts_by_nodenum = read_gmsh_nodes(io)

    # Read to the beginning of the elements section.
    if !read_through_line_starting(io, "\$Elements")
    error("Could not find beginning of elements section in mesh file.")
    end

    decl_num_els = uint64(readline(io)) # This count can include unwanted lower order elements.
    line, lines_read = read_until(io, is_polytope_or_endmarker)
    if line == "" || beginswith(line, "\$EndElements")
        error("No elements found in \$Elements section.")
    end
    est_base_tris = decl_num_els - (lines_read - 1)

    base_tris_iter = GmshElementsIter(line, io)

    (base_pts_by_nodenum, est_base_tris, base_tris_iter)

    faces = Array(Face,est_base_tris)

    i = 1
    for base_tri in base_tris_iter
        if i == length(faces)+1
            push!(faces,base_tri.point_nums)
        else
            faces[i] = base_tri.point_nums
        end
        i += 1
    end
    if i < length(faces)
        resize!(faces,i-1)
    end

    Mesh{Vertex3}(base_pts_by_nodenum,faces)
end

function read_through_line_starting(io::IO, line_start::ASCIIString)
    l = readline(io)
    while l != "" && !beginswith(l, line_start)
    l = readline(io)
    end
    l != ""
end

function read_gmsh_nodes(io::IO)
    if !read_through_line_starting(io, "\$Nodes")
        error("Nodes section not found in mesh input file.")
    else
        # Next line should be the vert count.
        const count = int(readline(io))
        const pts = Array(Vertex3, count)
        l = readline(io)
        while !beginswith(l, "\$EndNodes") && l != ""
            const toks = split(strip(l), ' ')
            if length(toks) != 4
                error("More than 4 columns for a node. Aborting!")
            end
            const pt = Vertex3(
                float64(toks[TOKN_NODELINE_POINT1]), 
                float64(toks[TOKN_NODELINE_POINT2]),
                float64(toks[TOKN_NODELINE_POINT3]))
            pts[uint64(toks[TOKN_NODELINE_ELNUM])] = pt
            l = readline(io)
        end
    end
    if !beginswith(l, "\$EndNodes") error("End of nodes section not found in mesh file.") end
    pts
end

function is_polytope_or_endmarker(line::ASCIIString)
    if beginswith(line, "\$EndElements")
        true
    else
        const toks = split(line, ' ')
        !is_lower_order_el_type(eltypenum(toks[TOKN_ELLINE_ELTYPE]))
    end
end

# Read the stream until the line condition function returns true or the stream is exhausted, returning
# either the line passing the condition, or "" if the stream was exhausted, and (in either case)
# the number of lines read including the successful line if any.
function read_until(io::IO, line_cond_fn::Function)
  l = readline(io)
  lines_read = 1
  while l != "" && !line_cond_fn(l)
    l = readline(io)
    lines_read += 1
  end
  l, lines_read
end

