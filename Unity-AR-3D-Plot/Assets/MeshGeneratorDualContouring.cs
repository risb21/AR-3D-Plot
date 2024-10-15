using System;
using System.Collections;
using System.Collections.Generic;
using System.Text.RegularExpressions;
using Unity.VisualScripting;
using UnityEngine;
using UnityEngine.Analytics;
using UnityEngine.Animations;
using UnityEngine.XR;

using cellIdx = System.ValueTuple<int, int, int>;
using coords = System.ValueTuple<float, float, float>;
using vector = System.ValueTuple<float, float, float>;

public delegate float SurfaceDelegate(coords p);

class DualContouring {

    // Stores all information needed for sampling data
    struct SampledData {
        public float resolution, bounds;
        public bool[,,] points;
    }

    struct Intersection {
        public coords point;
        public vector normal;
    };

    enum DC_Axis {
        X = 0,
        Y = 1,
        Z = 2,
    };

    SurfaceDelegate surface;
    SampledData samples;
    // public GameObject cubePrefab;

    static readonly cellIdx[] DC_Points = {
        (0, 0, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 0, 0),
        (0, 1, 0),
        (0, 1, 1),
        (1, 1, 1),
        (1, 1, 0),
    };
    static readonly (int, int)[] DC_Edges = {
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),

    };
    
    const int maxSamples = 128; // 2 ^ 8, i.e, 8 bit indexing
    const float delta = 1e-3f;
    
    public DualContouring(SurfaceDelegate surface, /* GameObject obj, */ float bounds = 5) {
        samples.bounds = bounds;
        // To include -ve & +ve bounds
        samples.resolution = 2 * bounds / (maxSamples - 1);
        this.surface = surface;
        // cubePrefab = obj;
    }

    void SamplePoints() {
        // Allocate memory for sampling
        samples.points = new bool[maxSamples, maxSamples, maxSamples];

        for (int i = 0; i < maxSamples; i++) {
            for (int j = 0; j < maxSamples; j++) {
                for (int k = 0; k < maxSamples; k++) {
                    // Convert indices to positions
                    float x = i * samples.resolution - samples.bounds;
                    float y = j * samples.resolution - samples.bounds;
                    float z = k * samples.resolution - samples.bounds;
                    
                    // Sample sign data
                    samples.points[i, j, k] = surface((x, y, z)) >= 0;
                }
            }
        }
    }

    coords UpdateIntersectionCoordinates(coords u, coords v, DC_Axis axis) {
        var (mx, my, mz) = u;
        var (ux, uy, uz) = u;
        var (vx, vy, vz) = v;

        switch (axis) {
            case DC_Axis.X:
                mx = (ux + vx) / 2;
                break;
            case DC_Axis.Y:
                my = (uy + vy) / 2;
                break;
            case DC_Axis.Z:
                mz = (uz + vz) / 2;
                break;
        }

        return (mx, my, mz);
    }

    DC_Axis GetAxisOfChange(cellIdx u, cellIdx v) {
        int axisChanges = 0;
        DC_Axis result = DC_Axis.X;
        
        var(ux, uy, uz) = u;
        var(vx, vy, vz) = v;

        if (ux != vx) {
            axisChanges++;
            result = DC_Axis.X;
        }

        if (uy != vy) {
            axisChanges++;
            result = DC_Axis.Y;
        }

        if (uz != vz) {
            axisChanges++;
            result = DC_Axis.Z;
        }

        Debug.Assert(axisChanges == 1);

        return result;
    }   

    // Takes 2 points on an edge, returns point of intersection of the surface on the edge
    coords ApproximateIntersection(cellIdx uIdx, cellIdx vIdx) {
        var (uIdxX, uIdxY, uIdxZ) = uIdx;
        var (vIdxX, vIdxY, vIdxZ) = vIdx;

        Debug.Assert(samples.points[uIdxX, uIdxY, uIdxZ] != samples.points[vIdxX, vIdxY, vIdxZ]);
        int maxIterations = 10;

        
        DC_Axis axis = GetAxisOfChange(uIdx, vIdx);

        // Ensure U indexes a negative sample
        if (samples.points[uIdxX, uIdxY, uIdxZ] == true) {
            switch (axis) {
                case DC_Axis.X:
                    (uIdxX, vIdxX) = (vIdxX, uIdxX);
                    break;
                case DC_Axis.Y:
                    (uIdxY, vIdxY) = (vIdxY, uIdxY);
                    break;
                case DC_Axis.Z:
                    (uIdxZ, vIdxZ) = (vIdxZ, uIdxZ);
                    break;
                
                // Impossible section
                default:
                    break;
            }
        }
        // // Detect axis of change in coords and Ensure U < V
        // if (uIdxX != vIdxX) {

        //     axis = DC_Axis.X;
        // } else if (uIdxY != vIdxY) {
        //     // Ensure U is smaller
        //     if (samples.points[uIdxX, uIdxY, uIdxZ] == false) {
        //         (uIdxY, vIdxY) = (vIdxY, uIdxY);
        //     }

        //     axis = DC_Axis.Y;
        // } else {
        //     // Ensure U is smaller
        //     if (samples.points[uIdxX, uIdxY, uIdxZ] == false) {
        //         (uIdxZ, vIdxZ) = (vIdxZ, uIdxZ);
        //     }

        //     axis = DC_Axis.Z;
        // }

        // Convert sample indices to coords 
        (float ux,float uy, float uz) = (
            uIdxX * samples.resolution - samples.bounds,
            uIdxY * samples.resolution - samples.bounds,
            uIdxZ * samples.resolution - samples.bounds
        );
        (float vx,float vy, float vz) = (
            vIdxX * samples.resolution - samples.bounds,
            vIdxY * samples.resolution - samples.bounds,
            vIdxZ * samples.resolution - samples.bounds
        );

        // Middle coords
        coords mid = (.0f, .0f, .0f);
        // Threshold below which value is considered zero
        float threshold = 1e-3f;
        
        // Binary Search for point of intersection on the edge
        while (maxIterations-- > 0) {

            // Update mid coords to be between U and V
            mid = UpdateIntersectionCoordinates((ux, uy, uz), (vx, vy, vz), axis);
            // Sample mid point
            float sampledPoint = surface(mid);

            if (Math.Abs(sampledPoint) < threshold) {
                break;
            }

            if (sampledPoint < 0) {
                // Move U ahead
                (ux, uy, uz) = mid;
            } else if (sampledPoint > 0) {
                // Move V behind
                (vx, vy, vz) = mid;
            }
        }

        return mid;
    }

    float GetGradient(coords pointOnSurface, DC_Axis axis) {
        var (x, y, z) = pointOnSurface;

        return axis switch {
            DC_Axis.X => (surface((x + delta, y, z)) - surface((x - delta, y, z))) / (2 * delta),
            DC_Axis.Y => (surface((x, y + delta, z)) - surface((x, y - delta, z))) / (2 * delta),
            DC_Axis.Z => (surface((x, y, z + delta)) - surface((x, y, z - delta))) / (2 * delta),
            // Impossible section
            _ => 0
        };
    }

    
    
    Dictionary<cellIdx, coords> CellTraverse() {

        int errors = 0;
        Dictionary<cellIdx, coords> cellCoords = new ();

        // i + 1 < maxSamples because no. of cells = no. of samples - 1
        for (int i = 0; i + 1 < maxSamples; i++) {
            for (int j = 0; j + 1 < maxSamples; j++) {
                for (int k = 0; k + 1 < maxSamples; k++) {
                    // 8 intersections per cell is a high amount
                    List<Intersection> intersections = new (8);

                    // The following will happen to every 3D cell
                    for (int edges = 0; edges < DC_Edges.Length; edges++) {
                        // Get points on edge
                        var (u, v) = DC_Edges[edges];
                        
                        // Get point offsets
                        var (ux, uy, uz) = DC_Points[u];
                        var (vx, vy, vz) = DC_Points[v];

                        // Get signs of points on the edge
                        bool uSign = samples.points[i + ux, j + uy, k + uz];
                        bool vSign = samples.points[i + vx, j + vy, k + vz];

                        if (uSign == vSign) { continue; }
                        
                        // Surface intersects edge when signs are different
                        coords surfacePt  = ApproximateIntersection(
                            (i + ux, j + uy, k + uz),
                            (i + vx, j + vy, k + vz)
                        );

                        // Find normal on the surface
                        // Slope on the surface, by displacing by tiny amount
                        coords surfaceNormal2 = (
                            GetGradient(surfacePt, DC_Axis.X),
                            GetGradient(surfacePt, DC_Axis.Y),
                            GetGradient(surfacePt, DC_Axis.Z)
                        );

                        Intersection currentCell;
                        currentCell.point = surfacePt;
                        currentCell.normal = surfaceNormal2;

                        intersections.Add(currentCell);
                    }

                    // Skip cell if there's zero intersections
                    if (intersections.Count == 0) {
                        continue;
                    }

                    // Bias the solve to Center of Mass between intersection points
                    // 1. Find CoM for points (average)
                    float comX = 0, comY = 0, comZ = 0;
                    foreach (Intersection intersection in intersections) {
                        var (pointX, pointY, pointZ) = intersection.point;
                        comX += pointX;
                        comY += pointY;
                        comZ += pointZ;
                    }

                    comX /= intersections.Count;
                    comY /= intersections.Count;
                    comZ /= intersections.Count;

                    // 2. Add bias vectors at CoM
                    const float bias = 0.25f;

                    Intersection xBias, yBias, zBias;
                    xBias.point = yBias.point = zBias.point = (comX, comY, comZ);
                    xBias.normal = (bias, .0f, .0f);
                    yBias.normal = (.0f, bias, .0f);
                    zBias.normal = (.0f, .0f, bias);
                    
                    intersections.Add(xBias);
                    intersections.Add(yBias);
                    intersections.Add(zBias);


                    // n normals, each of 3 values, (n x 3)
                    double[,] A = new double[intersections.Count, 3];
                    // bi = ni . pi , ni -> normal vector, pi -> point of intersection 
                    double[] b = new double[intersections.Count];

                    for (int l = 0; l < intersections.Count; l++) {
                        coords pl = intersections[l].point;
                        coords nl = intersections[l].normal;

                        var (plX, plY, plZ) = pl;
                        var (nlX, nlY, nlZ) = nl;

                        // Set normal matrix
                        A[l, (int) DC_Axis.X] = nlX;
                        A[l, (int) DC_Axis.Y] = nlY;
                        A[l, (int) DC_Axis.Z] = nlZ;

                        // Set bl = dot_product(nl, pl)
                        b[l] = plX * nlX + plY * nlY + plZ * nlZ;
                    }

                    // Solve for point that minimizes QEF
                    // Solve the system A*x = b
                    alglib.rmatrixsolvels(
                        A, intersections.Count, 3, b, 0.00000001,
                        out double[] x, out alglib.densesolverlsreport rep
                    );

                    if (rep.terminationtype != -4) {
                        cellCoords.Add(
                            (i, j, k),                                  // Cell index
                            ((float) x[0], (float) x[1], (float) x[2])  // Min QEF point
                        );

                        // Vector3 minQEF = new Vector3(
                        //     (float) x[0],
                        //     (float) x[2],
                        //     (float) x[1]
                        // );
                        
                        // Instantiate(cubePrefab, minQEF, Quaternion.identity);
                    } else {
                        errors += 1;
                    }
                }
            }
        }

        Debug.Log($"No. of solving errors: {errors}");

        return cellCoords;
    }

    (Vector3[], int[]) GenerateMesh(Dictionary<cellIdx, coords> cellsWithPoints) {
        List<Vector3> points = new ();
        List<int> tris = new ();

        // Mapping from Cell idx to Edge idx has to be consistent with all cells
        // Hence store edges as collections of 2 sample point indices
        // edge: (u, v), where sample(u) < 0 and u, v: coords
        (int xBound, int yBound, int zBound) = (
            samples.points.GetLength(0), 
            samples.points.GetLength(1),
            samples.points.GetLength(2) 
        );

        // Store edges which have been considered
        // NOTE: ONLY store with -ve sample first
        HashSet<(cellIdx, cellIdx)> completedEdges = new ();

        foreach (KeyValuePair<cellIdx, coords> entry in cellsWithPoints) {
            var (cellX, cellY, cellZ) = entry.Key;

            // Consider all intersections on edges for each cell with a point
            for (int edge = 0; edge < DC_Edges.Length; edge++) {
                // Get point indices on the edge
                var (u, v) = DC_Edges[edge];

                // Get sample point offsets wrt cell idx
                var (ux, uy, uz) = DC_Points[u];
                var (vx, vy, vz) = DC_Points[v];

                var (uPtX, uPtY, uPtZ) = (cellX + ux, cellY + uy, cellZ + uz);
                var (vPtX, vPtY, vPtZ) = (cellX + vx, cellY + vy, cellZ + vz);

                var uSign = samples.points[uPtX, uPtY, uPtZ];
                var vSign = samples.points[vPtX, vPtY, vPtZ];

                // Check if edge has intersection
                if (uSign == vSign) { continue; }

                // Ensure U is the smaller index
                if (uPtX > vPtX || uPtY > vPtY || uPtZ > vPtZ) {
                    (uPtX, uPtY, uPtZ, vPtX, vPtY, vPtZ) = (vPtX, vPtY, vPtZ, uPtX, uPtY, uPtZ);
                    (uSign, vSign) = (vSign, uSign);
                }

                // Check if edge has been considered already
                if (completedEdges.Contains(((uPtX, uPtY, uPtZ), (vPtX, vPtY, vPtZ)))) {
                    continue;
                }

                completedEdges.Add(((uPtX, uPtY, uPtZ), (vPtX, vPtY, vPtZ)));

                // Determine axis of change along edge
                DC_Axis changeAxis = GetAxisOfChange((uPtX, uPtY, uPtZ), (vPtX, vPtY, vPtZ));

                // Ensure coords are within bounds on  
                // axes other than the axis of change
                bool withinBounds = true;
                foreach (DC_Axis axis in Enum.GetValues(typeof(DC_Axis))) {
                    if (axis == changeAxis) { continue; }

                    switch (axis) {
                        case DC_Axis.X:
                            // Only checking one point because both 
                            // points must have same indices at axes 
                            // other than the axis of change
                            if (uPtX == 0 || uPtX == xBound - 1) {
                                withinBounds = false;
                            }
                            break;
                        case DC_Axis.Y:
                            if (uPtY == 0 || uPtY == yBound - 1) {
                                withinBounds = false;
                            }
                            break;
                        case DC_Axis.Z:
                            if (uPtZ == 0 || uPtZ == zBound - 1) {
                                withinBounds = false;
                            }
                            break;
                    }
                }

                if (!withinBounds) { continue; }
                
                // Gather points around the intersected edge

                // U is the smaller index, so all 4 points 
                // can be gathered by offsetting it
                cellIdx[] offsets = new cellIdx[4];
                offsets[0] = ( 0,  0,  0 );

                // These are clockwise if sign at U is +ve
                switch (changeAxis) {
                    case DC_Axis.X:
                        offsets[1] = ( 0, -1,  0 );
                        offsets[2] = ( 0, -1, -1 );
                        offsets[3] = ( 0,  0, -1 );
                        break;
                    case DC_Axis.Y:
                        offsets[1] = ( 0,  0, -1 );
                        offsets[2] = (-1,  0, -1 );
                        offsets[3] = (-1,  0,  0 );
                        break;
                    case DC_Axis.Z:
                        offsets[1] = (-1,  0,  0 );
                        offsets[2] = (-1, -1,  0 );
                        offsets[3] = ( 0, -1,  0 );
                        break;
                }
                
                List<Vector3> tempPoints = new (4);
                
                foreach (cellIdx offset in offsets) {
                    var (offsetX, offsetY, offsetZ) = offset;

                    if (!cellsWithPoints.ContainsKey((
                        uPtX + offsetX,
                        uPtY + offsetY,
                        uPtZ + offsetZ
                    ))) { continue; }

                    var (vertX, vertY, vertZ) = cellsWithPoints[(
                        uPtX + offsetX,
                        uPtY + offsetY,
                        uPtZ + offsetZ
                    )];

                    tempPoints.Add(new Vector3(vertX, vertZ, vertY));
                }

                if (tempPoints.Count <= 2) { continue; }

                // Determine winding direction
                if (!uSign) {
                    tempPoints.Reverse();
                }

                // Adding 1 or 2 tris (depending on no. of points found)
                // Outside faces first
                points.Add(tempPoints[0]); tris.Add(points.Count - 1);
                points.Add(tempPoints[1]); tris.Add(points.Count - 1);
                points.Add(tempPoints[2]); tris.Add(points.Count - 1);
                
                if (tempPoints.Count == 4) {
                    points.Add(tempPoints[2]); tris.Add(points.Count - 1);
                    points.Add(tempPoints[3]); tris.Add(points.Count - 1);
                    points.Add(tempPoints[0]); tris.Add(points.Count - 1);
                }

                // Now inside faces 
                points.Add(tempPoints[2]); tris.Add(points.Count - 1);
                points.Add(tempPoints[1]); tris.Add(points.Count - 1);
                points.Add(tempPoints[0]); tris.Add(points.Count - 1);
                
                if (tempPoints.Count == 4) {
                    points.Add(tempPoints[0]); tris.Add(points.Count - 1);
                    points.Add(tempPoints[3]); tris.Add(points.Count - 1);
                    points.Add(tempPoints[2]); tris.Add(points.Count - 1);
                }

            }
        }

        return (points.ToArray(), tris.ToArray());
    }

    public (Vector3[], int[]) Run() {
        SamplePoints();
        var cellsToPoints = CellTraverse();
        
        return GenerateMesh(cellsToPoints);
    }
}


[RequireComponent(typeof(MeshFilter))]
public class MeshGeneratorDualContouring : MonoBehaviour {

    Mesh mesh;
    // public GameObject pointPrefab;
    Vector3[] verts;
    int[] tris;

    int r = 1;
    readonly float R = 2f;


    void Awake() {
        mesh = GetComponent<MeshFilter>().mesh;
    }

    float Surface(coords p) {
        var (x, y, z) = p;
        
        // Torus
        // return z * z + (float) Math.Pow(R - Math.Sqrt(x * x + y * y), 2) - r * r; 
        
        // Sphere
        // return z * z + y * y + x * x - R * R;

        // Implicit surface
        // return x * x + y * y + z * z + (float) (Math.Sin(4 * x) + Math.Sin(4 * y) + Math.Sin(4 * z)) - R;
        
        // Idk
        return 4 * (float) Math.Pow(Math.E, -y * y / 4) * (float) Math.Sin(2 * x) - z;
        
        // Goursat's tangle
        // return (float) (Math.Pow(x, 4) + Math.Pow(y, 4) + Math.Pow(z, 4)) / 2
        //        - (x * x + y * y + z * z) * 8 + 60;

        // With a cusp
        // return (float) (Math.Pow(x, 2/3) + Math.Pow(y, 2/3)) - z;
       
    }

    void Start() {
        // GenerateMeshData();
        // CreateMesh();

        DualContouring d = new DualContouring(Surface, 6f);

        DateTime start = DateTime.Now;

        mesh.Clear();
        mesh.indexFormat = UnityEngine.Rendering.IndexFormat.UInt32;
        (mesh.vertices, mesh.triangles) = d.Run();

        mesh.RecalculateNormals();

        Debug.Log($"Time taken: {DateTime.Now - start}");
    }

    void GenerateMeshData() {
        verts = new Vector3[]{new Vector3(0, 0, 0), new Vector3(0, 0, 1), new Vector3(1, 0, 0)};
        tris = new int[]{0, 1, 2};
    }
    
    void CreateMesh() {
        mesh.Clear();
        mesh.vertices = verts;
        mesh.triangles = tris;
        mesh.RecalculateNormals();
    }

}
