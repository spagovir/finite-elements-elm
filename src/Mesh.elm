module Mesh exposing (..)
import Array
import Maybe
import List

sequenceMaybe : List (Maybe.Maybe a) -> Maybe.Maybe (List a)
sequenceMaybe = List.foldr (Maybe.map2 (::)) (Just [])


type alias MeshPoint = 
        { p : Vec2
        , angles : List (Int, Int) -- meshpoints are referenced by their array index. 
        }
-- Meshes are stored as adjacency arrays, except
-- instead of storing edges we store angles (pairs of adjacent edges)
type alias Mesh = Array.Array MeshPoint 


type alias Vec2 =
       { x : Float
       , y : Float
       } 

add : Vec2 -> Vec2 -> Vec2
add v1 v2 = Vec2 (v1.x + v2.x) (v1.y + v2.y)

sub : Vec2 -> Vec2 -> Vec2
sub v1 v2 = Vec2 (v1.x - v2.x) (v1.y - v2.y)

dot : Vec2 -> Vec2 -> Float
dot v1 v2 = v1.x * v2.x + v1.y * v2.y 

scale : Float -> Vec2 -> Vec2
scale k v1 = Vec2 (k * v1.x) (k* v1.y)

cross : Vec2 -> Vec2 -> Float
cross v1 v2 = v1.x * v2.y - v1.y * v2.x


-- Gives the map from R^n to itself through H_0 then H_0^* then (R^n)^* by
-- identifying a vector x in R^n with the values of a piecwise linear
-- function u_x such that if $p_i$ is the meshpoint with index i,
-- u_x(p_i) = u_i.
--
-- That is to say, we are using the mesh to create the linear transform T 
-- the property that, $T(x) . y = int_omega u_x * u_y$ for all $y \in R^n$.  
--
-- The function takes an additional integer parameter m specifying that the
-- first m meshpoints are boundary points and are therefore constrained to be
-- 0 under the Dirichlet boundary condition.
--
-- Both the input and output vectors should have the same length as the number 
-- of meshpoints in the mesh. 
--
-- Fails if graph is ill defined (ie, contains pointers to vertices that do not exist)
-- or if # of vertices in graph and dimension of vector does not match. 
toFnl : Mesh -> Int -> Array.Array Float -> Maybe (Array.Array Float)
toFnl mesh m x = 
        List.range 0 (Array.length mesh - 1)
                |> List.map (\k -> if k < m then Just 0 else testElem mesh x k)
                |> sequenceMaybe
                |> Maybe.map Array.fromList
triangInnerProduct mesh x k0 ks =
        let (k1,k2) = ks in
        (Array.get k0 x)
               |> Maybe.andThen (\x0 -> Array.get k1 x
               |> Maybe.andThen (\x1 -> Array.get k2 x
               |> Maybe.andThen (\x2 -> Array.get k0 mesh
               |> Maybe.andThen (\p0 -> Array.get k1 mesh
               |> Maybe.andThen (\p1 -> Array.get k2 mesh
               |> Maybe.map (\p2 -> 
                       let area4 = 2.0 * abs ( cross (sub p1.p p2.p) (sub p2.p p0.p)) in
                           x0 * dot (sub p1.p p2.p) (sub p1.p p2.p) / area4
                                + x1 * dot (sub p2.p p0.p) (sub p2.p p1.p) / area4 
                                + x2 * dot (sub p1.p p0.p) (sub p1.p p2.p) / area4))))))
testElem mesh x k =
               (Array.get k mesh)
                      |> Maybe.map .angles
                      |> Maybe.map (List.map (triangInnerProduct mesh x k))
                      |> Maybe.andThen sequenceMaybe
                      |> Maybe.map (List.sum)
                    
-- Pointer used to index a member of a disjoint union of two arrays.
type DisjointPtr = B Int | I Int 

type alias Triangle = {edge1 : DisjointPtr, edge2 : DisjointPtr, edge3 : DisjointPtr}

-- An intermediate representation of a triangular mesh for computing mesh refinements.
type alias IntermediateMesh = 
        { bPoints : Int
        , iPoints : Int
        , bEdges : Array.Array (DisjointPtr, DisjointPtr) -- ptrs are assumed to be in order, with all boundary vertices < internal vertices
        , iEdges : Array.Array (DisjointPtr, DisjointPtr)
        , triangles : Array.Array Triangle } -- triangles are three edges in lexicographic order

-- Finds the vertex pointer in the new mesh corresponding to the midpoint of an edge in the old mesh. 
midpoint : IntermediateMesh -> DisjointPtr -> DisjointPtr
midpoint mesh edge =
        case edge of 
                B k -> B (mesh.bPoints + k)
                I k -> I (mesh.iPoints + k)



refineBPoints mesh = mesh.bPoints + Array.length mesh.bEdges
refineIPoints mesh = mesh.iPoints + Array.length mesh.iEdges
refineBEdges mesh = 
        Array.toIndexedList mesh.bEdges
        |> List.concatMap (\(k, edge) -> [(Tuple.first edge, midpoint mesh (B k)),(Tuple.second edge, midpoint mesh (B k))])
        |> Array.fromList

refineIEdges mesh =
        let
            fromEdges = 
                    Array.toIndexedList mesh.iEdges 
                    |> List.concatMap (\(k, edge) -> [(Tuple.first edge, midpoint mesh (I k)), (Tuple.second edge, midpoint mesh (I k))]) 
            fromTriangs =
                    Array.toIndexedList mesh.triangles
                    |> List.concatMap (\(k, triang) -> [
                                    (midpoint mesh (triang.edge1), midpoint mesh (triang.edge2)),
                                    (midpoint mesh (triang.edge1), midpoint mesh (triang.edge3)),
                                    (midpoint mesh (triang.edge2), midpoint mesh (triang.edge3))])
        in  
            Array.fromList (fromEdges ++ fromTriangs)

subEdgePtrs : DisjointPtr -> (DisjointPtr, DisjointPtr)
subEdgePtrs edge =
        case edge of
                B k -> (B (2 * k), B (2 * k + 1))
                I k -> (I (2 * k), I (2 * k + 1))

triangEdgePtrs : IntermediateMesh -> Int -> Triangle 
triangEdgePtrs mesh k = let t = 2 * Array.length mesh.iEdges + 3 * k in Triangle (I t) (I (t + 1)) (I (t + 2))

refineTriangs mesh =
        Array.toIndexedList mesh.triangles
        |> List.concatMap (\(k, triang) -> [
                Triangle (Tuple.first (subEdgePtrs triang.edge1)) (Tuple.first (subEdgePtrs triang.edge2)) ((triangEdgePtrs mesh k).edge1),
                Triangle (Tuple.second (subEdgePtrs triang.edge1)) (Tuple.first (subEdgePtrs triang.edge3)) ((triangEdgePtrs mesh k).edge2),
                Triangle (Tuple.second (subEdgePtrs triang.edge2)) (Tuple.second (subEdgePtrs triang.edge3)) ((triangEdgePtrs mesh k).edge3),
                triangEdgePtrs mesh k])
        |> Array.fromList

refine : IntermediateMesh -> IntermediateMesh
refine mesh = IntermediateMesh (refineBPoints mesh) (refineIPoints mesh) (refineBEdges mesh) (refineIEdges mesh) (refineTriangs mesh)



