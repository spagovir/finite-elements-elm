module Main exposing (..)
import Array
import Array.Extra
import Maybe
import List
import Svg as S
import Svg.Attributes as SA
import Hex
import Result
import Browser
import Html exposing (..)

sequenceMaybe : Array.Array (Maybe.Maybe a) -> Maybe.Maybe (Array.Array a)
sequenceMaybe = Array.foldl (Maybe.map2 (Array.push)) (Just Array.empty)

sequenceMaybeL : List (Maybe.Maybe a) -> Maybe.Maybe (List a)
sequenceMaybeL = List.foldr (Maybe.map2 (::)) (Just [])


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
        Array.initialize (Array.length x) (\k -> if k < m then Just 0 else testElem mesh x k)
        |> sequenceMaybe
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
                                - x1 * dot (sub p2.p p0.p) (sub p2.p p1.p) / area4 
                                - x2 * dot (sub p1.p p0.p) (sub p1.p p2.p) / area4))))))
testElem mesh x k =
               (Array.get k mesh)
                      |> Maybe.map .angles
                      |> Maybe.map (List.map (triangInnerProduct mesh x k))
                      |> Maybe.andThen sequenceMaybeL
                      |> Maybe.map (List.sum)

aDot : Array.Array Float -> Array.Array Float -> Float
aDot a b = Array.Extra.map2 (*) a b |> (Array.foldr (+) 0)

aSub : Array.Array Float -> Array.Array Float -> Array.Array Float
aSub = Array.Extra.map2 (-)

aScale : Float -> Array.Array Float -> Array.Array Float
aScale = Array.map << (*)

cgStep : (Array.Array Float -> Maybe (Array.Array Float)) -> Array.Array Float -> Array.Array Float -> Array.Array Float -> Array.Array Float -> Int -> Maybe (Array.Array Float)
cgStep t u b r s n =
        if n == 0 
        then Just u
        else
                t s 
                |> Maybe.andThen (\ts -> aSub u (aScale ((aDot r s) / (aDot s ts)) s)
                |> (\nu -> t nu
                |> Maybe.andThen (\tu -> aSub tu b
                |> (\nr -> cgStep t nu b nr (aSub nr (aScale (aDot nr ts / aDot s ts) s)) (n-1)))))

cg : (Array.Array Float -> Maybe (Array.Array Float)) -> Array.Array Float -> Array.Array Float -> Int -> Maybe (Array.Array Float)
cg t u b n =
        t u
        |> Maybe.andThen (\tu -> cgStep t u b (aSub tu b) (aSub tu b) n)

gradDescent t u b n =
        if n == 0
        then Just u
        else
                t u
                |> Maybe.andThen (\tu -> aSub tu b
                |> (\r -> t r
                |> Maybe.andThen (\tr -> gradDescent t (aSub u (aScale (aDot r r / aDot tr r) r)) b (n-1)))) 
                    

testMatrix v = Array.fromList [aDot (Array.fromList [1, 2]) v, aDot (Array.fromList [3,4]) v] |> Just
        -- Array.fromList [aDot (Array.fromList [1,2,3,4]) v, aDot (Array.fromList [5,6,7,8]) v, aDot (Array.fromList [9,10,11,12]) v, aDot (Array.fromList [13,14,15,16]) v] |> Just 

laplaceVector : Mesh -> Int -> Float -> Maybe (Array.Array Float)
laplaceVector mesh m l = 
        Array.initialize 
        (Array.length mesh) 
        (\k -> 
                if k < m 
                then Just 0.0 
                else 
                        Array.get k mesh 
                        |> Maybe.andThen (\mpoint -> 
                                mpoint.angles 
                                |> List.map (area mesh mpoint) 
                                |> sequenceMaybeL 
                                |> Maybe.map (List.sum >> (*) (-l/3) ))) 
        |> sequenceMaybe 

area : Mesh -> MeshPoint -> (Int, Int) -> Maybe Float
area mesh p0 ps = 
        Array.get (Tuple.first ps) mesh 
        |> Maybe.andThen (\p1 -> Array.get (Tuple.second ps) mesh
        |> Maybe.map (\p2 -> abs ( cross (sub p1.p p2.p) (sub p2.p p0.p)) / 2.0)) 
              
       
                    
-- Pointer used to index a member of a disjoint union of two arrays.
type DisjointPtr = B Int | I Int 

type alias Triangle = {edge1 : DisjointPtr, edge2 : DisjointPtr, edge3 : DisjointPtr}

-- An intermediate representation of a triangular mesh for computing mesh refinements.
type alias IntermediateMesh = 
        { bPoints : Array.Array Vec2
        , iPoints : Array.Array Vec2
        , bEdges : Array.Array (DisjointPtr, DisjointPtr) -- ptrs are assumed to be in order, with all boundary vertices < internal vertices
        , iEdges : Array.Array (DisjointPtr, DisjointPtr)
        , triangles : Array.Array Triangle } -- triangles are three edges in lexicographic order

disjointGet : DisjointPtr -> IntermediateMesh -> Maybe Vec2 -- gets a point referenced by a disjoint ptr. 
disjointGet i mesh =
        case i of 
                B k -> Array.get k mesh.bPoints
                I k -> Array.get k mesh.iPoints

disjointGetEdges i mesh = 
        case i of 
                B k -> Array.get k mesh.bEdges
                I k -> Array.get k mesh.iEdges
-- gets a point referenced by a disjoint ptr from an array that has the length of bPoints ++ iPoints
concatGet : DisjointPtr -> IntermediateMesh -> Array.Array a -> Maybe a
concatGet i mesh array =
        case i of
                B k -> Array.get k array
                I k -> Array.get (k + Array.length mesh.bPoints) array

-- Finds the vertex pointer in the new mesh corresponding to the midpoint of an edge in the old mesh. 
midpoint : IntermediateMesh -> DisjointPtr -> DisjointPtr
midpoint mesh edge =
        case edge of 
                B k -> B (Array.length mesh.bPoints + k)
                I k -> I (Array.length mesh.iPoints + k)



refineBPoints mesh = 
        mesh.bEdges
        |> Array.map (\(ptr1,ptr2) ->
                disjointGet ptr1 mesh 
                |> Maybe.andThen (\v1 -> disjointGet ptr2 mesh
                |> Maybe.map (\v2 -> add v1 v2 |> scale 0.5))) 
        |> sequenceMaybe
        |> Maybe.map (Array.append mesh.bPoints)
refineIPoints mesh =
        mesh.iEdges
        |> Array.map (\(ptr1,ptr2) ->
                disjointGet ptr1 mesh
                |> Maybe.andThen (\v1 -> disjointGet ptr2 mesh
                |> Maybe.map (\v2 -> add v1 v2 |> scale 0.5))) 
        |> sequenceMaybe
        |> Maybe.map (Array.append mesh.iPoints)
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

refine : IntermediateMesh -> Maybe IntermediateMesh
refine mesh = 
        refineBPoints mesh
        |> Maybe.andThen (\nbPoints -> refineIPoints mesh
        |> Maybe.map (\niPoints -> IntermediateMesh nbPoints niPoints (refineBEdges mesh) (refineIEdges mesh) (refineTriangs mesh)))

-- Given a mesh of triangles and an array of temperatures at meshpoints, create a heat map 
-- with gradient between the temperatures of the corners on every triangle. 
coordString : Vec2 -> String 
coordString v = 
        String.fromFloat v.x ++ "," ++ String.fromFloat v.y

drawTriangleGradient : IntermediateMesh -> Array.Array Float -> (Float -> String) -> (Int, Triangle) -> Maybe (List (S.Svg MSG))
drawTriangleGradient mesh temps colorMap ktriang =
        let (k, triang) = ktriang in 
        disjointGetEdges triang.edge1 mesh
        |> Maybe.andThen (\(p1,p2) -> disjointGetEdges triang.edge2 mesh
        |> Maybe.andThen (\(_,p3) -> disjointGet p1 mesh
        |> Maybe.andThen (\p1x -> disjointGet p2 mesh
        |> Maybe.andThen (\p2x -> disjointGet p3 mesh
        |> Maybe.andThen (\p3x -> concatGet p1 mesh temps 
        |> Maybe.andThen (\p1t -> concatGet p2 mesh temps
        |> Maybe.andThen (\p2t -> concatGet p3 mesh temps
        |> Maybe.map (\p3t -> 
                let ((q1t,q1x),(q2t,q2x),(q3t,q3x)) = sortBy3 Tuple.first ((p1t,p1x),(p2t,p2x),(p3t,p3x))
                    gradBase = 
                            if q2t < q1t 
                            then 
                                    let 
                                        q4x = scale ((q1t-q3t) / (q1t - q2t)) (sub q2x q1x) |> add q1x  
                                        s = dot (sub q1x q3x) (sub q4x q3x) / dot (sub q4x q3x) (sub q4x q3x) 
                                    in 
                                    scale s (sub q4x q3x) |> add q3x 
                            else if q3t < q1t then sub q3x (scale (dot (sub q3x q1x) (sub q2x q1x) / dot (sub q2x q1x) (sub q2x q1x)) (sub q2x q1x)) 
                            else q1x

                               
                    linGrad =
                            S.linearGradient
                            [ SA.id ("grad" ++ String.fromInt k)
                            , SA.x1 (String.fromFloat q1x.x)
                            , SA.y1 (String.fromFloat q1x.y)
                            , SA.x2 (String.fromFloat gradBase.x)
                            , SA.y2 (String.fromFloat gradBase.y)
                            , SA.gradientUnits "userSpaceOnUse"
                            ]
                            [ S.stop [SA.stopColor (colorMap q1t), SA.offset "0%"] []
                            , S.stop [SA.stopColor (colorMap q3t), SA.offset "100%"] []
                            ]
                    triangSVG = 
                            S.polygon
                            [ SA.points ((coordString p1x) ++ " " ++ coordString p2x ++ " " ++ coordString p3x)
                            , SA.fill ("url('#grad" ++ String.fromInt k ++ "')")
                            ]
                            []
                in [linGrad, triangSVG]))))))))
                
sortBy3 : (a -> comparable) -> (a,a,a) -> (a,a,a) -- sorts three elements assending
sortBy3 f is =
        let (i1, i2, i3) = is
        in 
        if (f i1 >= f i2 && f i2 >= f i3) then (i1, i2, i3)
        else if (f i1 >= f i3 && f i3 >= f i2) then (i1, i3, i2)
        else if (f i2 >= f i1 && f i1 >= f i3) then (i2, i1, i3)
        else if (f i2 >= f i3 && f i3 >= f i1) then (i2, i3, i1)
        else if (f i1 >= f i2) then (i3 , i1, i2)
        else (i3, i2, i1)


        

heatMap : IntermediateMesh -> Array.Array Float -> Maybe (S.Svg MSG)
heatMap mesh temps =
        let 
            -- we want to map the hottest meshpoint to yellow and the coldest to blue. 
            colorMap t =
                   let maxTemp = Maybe.withDefault 1.0 (List.maximum (Array.toList temps))
                       minTemp = Maybe.withDefault 0.0 (List.minimum (Array.toList temps))
                       r = Hex.toString (floor ((t-minTemp)/(maxTemp-minTemp) * 255.5)) |> String.padLeft 2 '0'
                       g = Hex.toString (floor ((t-minTemp)/(maxTemp-minTemp) * 255.5)) |> String.padLeft 2 '0'
                       b = Hex.toString (floor ((maxTemp - t)/(maxTemp-minTemp) * 128.0)) |> String.padLeft 2 '0'
                   in "#" ++ r ++ g ++ b
        in
        mesh.triangles
        |> Array.toIndexedList 
        |> List.map (drawTriangleGradient mesh temps colorMap)
        |> sequenceMaybeL
        |> Maybe.map List.concat
        |> Maybe.map (S.svg [ SA.width "800", SA.height "800", SA.viewBox "0 0 800 800"])

toRegPtr mesh ptr = 
        case ptr of 
                B k -> k
                I k -> k + Array.length (mesh.bPoints)

addTriangle : IntermediateMesh -> Triangle -> Mesh -> Maybe Mesh 
addTriangle imesh triang mesh =
        disjointGetEdges triang.edge1 imesh
        |> Maybe.andThen (\(p1,p2) -> disjointGetEdges triang.edge2 imesh
        |> Maybe.andThen (\(_,p3) -> addAngles imesh mesh p1 p2 p3)) 

mutate : Int -> (a -> a) -> Array.Array a -> Maybe (Array.Array a)
mutate i f arr = Array.get i arr |> Maybe.map (\x -> Array.set i (f x) arr)

addAngles imesh mesh p1 p2 p3 =
      let 
          p1r = toRegPtr imesh p1  
          p2r = toRegPtr imesh p2
          p3r = toRegPtr imesh p3
      in
      mutate p1r (\mp -> MeshPoint (mp.p) ((p2r,p3r) :: mp.angles)) mesh 
      |> Maybe.andThen (\mesh1 -> mutate p2r (\mp -> MeshPoint (mp.p) ((p1r,p3r) :: mp.angles)) mesh1
      |> Maybe.andThen (\mesh2 -> mutate p3r (\mp -> MeshPoint (mp.p) ((p1r,p2r) :: mp.angles)) mesh2 ))

toAngleMesh : IntermediateMesh -> Maybe Mesh
toAngleMesh mesh = Array.foldr (Maybe.andThen << addTriangle mesh) (Just (Array.map (\v -> MeshPoint v []) (Array.append mesh.bPoints mesh.iPoints))) mesh.triangles 

calcTemps : Float -> IntermediateMesh -> Maybe (Array.Array Float)
calcTemps l imesh = 
        toAngleMesh imesh 
        |> Maybe.andThen (\mesh -> laplaceVector mesh (Array.length imesh.bPoints) l
        |> Maybe.andThen (\lvec -> cg (toFnl mesh (Array.length imesh.bPoints)) (Array.repeat (Array.length mesh) 0.0) lvec (1*(Array.length mesh - (Array.length imesh.bPoints)))))


main = -- gradDescent testMatrix (Array.fromList [0, 0]) (Array.fromList [3, 7]) 20 |> Maybe.withDefault (Array.fromList [0]) |> Array.map (String.fromFloat) |> Array.toList |> String.join "," |> text 
       -- cg testMatrix (Array.repeat 4 0) (Array.fromList [10,26,42,58]) 4 |> Maybe.withDefault (Array.fromList [0,0,0,0]) |> Array.map (String.fromFloat) |> Array.toList |> String.join ", " |> text 
       Browser.element {init=init, update=update, subscriptions = always Sub.none, view=view}

type MSG = 
        ToggleMode
        | Point DisjointPtr
        | Edge ((DisjointPtr, DisjointPtr), DisjointPtr)
        | Background Vec2

edgeToRegPtrs : IntermediateMesh -> (DisjointPtr,DisjointPtr) -> (Int, Int)
edgeToRegPtrs mesh = Tuple.mapBoth (toRegPtr mesh) (toRegPtr mesh)

type alias Model = 
        { shapes : Maybe IntermediateMesh
        , mode : Mode
        }

type Mode = Draw Tool | Display (Maybe IntermediateMesh) (Maybe (Array.Array Float))
type Tool = NoTool | BPoint | IPoint | BEdge1 | BEdge2 DisjointPtr | IEdge1 | IEdge2 DisjointPtr | Triang1 | Triang2 ((DisjointPtr, DisjointPtr),DisjointPtr) | Triang3 ((DisjointPtr,DisjointPtr), DisjointPtr) ((DisjointPtr,DisjointPtr),DisjointPtr)

drawWithError : Result.Result String (Html MSG) -> Html MSG
drawWithError result = 
        case result of 
                Ok okResult -> okResult
                Err errMsg -> text errMsg

view : Model -> Html MSG
view model =
        case model.mode of 
                Display shapes temps -> 
                        Result.fromMaybe "Error building mesh -> likely bad graph pointer" shapes 
                        |> Result.andThen (\jShapes -> Result.fromMaybe "Error calculating temperatures -> bad graph or vector length mismatch?" temps
                        |> Result.andThen (\jTemps -> Result.fromMaybe "Error matching temperatures to meshpoints" (heatMap jShapes jTemps)))
                        |> drawWithError
                Draw tool -> text "Draw mode under construction."

update : MSG -> Model -> (Model, Cmd MSG)
update msg model = (model, Cmd.none)

calcDrawMode : IntermediateMesh -> Mode
calcDrawMode imesh =
        let refinedMesh = refineN 4 imesh
        in
        Display refinedMesh (refinedMesh |> Maybe.andThen (calcTemps -1.0))

refineN : Int -> IntermediateMesh -> Maybe IntermediateMesh
refineN n mesh = 
        if n == 0 then Just mesh else refine mesh |> Maybe.andThen (refineN (n-1)) 

init : () -> (Model, Cmd MSG)
init _ = 
       let
          simpleTriang = 
                  IntermediateMesh 
                  ( Array.fromList [Vec2 0.0 0.0, Vec2 800.0 0.0, Vec2 0.0 800.0])
                  ( Array.empty)
                  ( Array.fromList [(B 0, B 1), (B 0, B 2), (B 1, B 2)])
                  ( Array.empty)
                  ( Array.fromList [Triangle (B 0) (B 1) (B 2)] )
          incisedSquare =
                  IntermediateMesh
                  ( Array.fromList [Vec2 0.0 0.0, Vec2 600.0 0.0, Vec2 600.0 600.0, Vec2 0.0 600.0, Vec2 300.0 300.0])
                  ( Array.empty)
                  ( Array.fromList [(B 0, B 1), (B 0, B 3), (B 0, B 4), (B 1, B 2), (B 2, B 3)])
                  ( Array.fromList [(B 1, B 4), (B 2, B 4), (B 3, B 4)])
                  ( Array.fromList [Triangle (B 0) (B 2) (I 0), Triangle (B 1) (B 2) (I 2), Triangle (B 3) (I 0) (I 1), Triangle (B 4) (I 1) (I 2)])
       in
       (Model (Just incisedSquare) (calcDrawMode incisedSquare), Cmd.none)
                    



