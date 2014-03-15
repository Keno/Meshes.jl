using ImmutableArrays

abstract Vertex

import Base: getindex

immutable Vertex3 <: Vertex
    coords::Vector3{Float64}
    Vertex3(a,b,c) = new(Vector3{Float64}(a,b,c))
    Vertex3(a::Vector3{Float64}) = new(a)
end

immutable Vertex2 <: Vertex
    coords::Vector2{Float64}
    Vertex2(a,b) = new(Vector2{Float64}(a,b))
    Vertex2(a::Vector2{Float64}) = new(a)
end

getindex(v::Vertex,args...) = getindex(v.coords,args...)

immutable Face
    v1 :: Int64
    v2 :: Int64
    v3 :: Int64
end

abstract AbstractMesh

type Mesh{VertexT} <: AbstractMesh
    vertices :: Vector{VertexT}
    faces :: Vector{Face}
end
#Mesh{VertexT}(a::Vector{VertexT},b::Vector{Face}) = Mesh{VertexT}(a,b)

vertices(m::Mesh) = m.vertices
faces(m::Mesh) = m.faces

# concatenates two meshes
function merge(m1::AbstractMesh, m2::AbstractMesh)
    v1 = vertices(m1)
    f1 = faces(m1)
    v2 = vertices(m2)
    f2 = faces(m2)
    nV = size(v1,1)
    nF = size(f2,1)
    newF2 = Face[ Face(f2[i].v1+nV, f2[i].v2+nV, f2[i].v3+nV) for i = 1:nF ]
    Mesh(append!(v1,v2),append!(f1,newF2))
end

export Vertex2, Vertex3, Face, AbstrctMesh, Mesh, vertices, faces, merge
