package main

import (
	"fmt"
	"image"
	"image/color"
	"image/draw"
	"image/png"
	"io/ioutil"
	"log"
	"math"
	"os"
	"time"

	"github.com/go-gl/mathgl/mgl32"
	"github.com/udhos/gwob"
)

// MaxDist tells us how far to check for things
const MaxDist = 1000

// ImgWidth width of the output
const ImgWidth = 256

// ImgHeight height of the output
const ImgHeight = 256

// Ray describes a ray starting from a position and moving in a direction
type Ray struct {
	Origin    mgl32.Vec3
	Direction mgl32.Vec3
}

// Hit provides information on where a ray hits something
type Hit struct {
	Length, Index float32
}

func (r *Ray) String() string {
	return fmt.Sprintf("%v to %v", r.Origin, r.Direction)
}

// Triangle is defined by 3 points
type Triangle struct {
	A, B, C mgl32.Vec3
}

// Normal determines the triangle's face normal
func (t *Triangle) Normal() mgl32.Vec3 {
	ab := t.B.Sub(t.A)
	ac := t.C.Sub(t.A)
	un := ab.Cross(ac)
	return un.Normalize()
}

// IntersectRay determines whether a ray hits the triangle. moller - trumbore
func (t *Triangle) IntersectRay(r *Ray) (bool, Hit) {
	ab := t.B.Sub(t.A)
	ac := t.C.Sub(t.A)
	hh := r.Direction.Cross(ac)
	aa := ab.Dot(hh)
	eps := float32(0.0001)
	if aa > -eps && aa < eps {
		return false, Hit{}
	}
	ff := 1 / aa
	ss := r.Origin.Sub(t.A)
	uu := ss.Dot(hh) * ff
	if uu < 0 || uu > 1 {
		return false, Hit{}
	}
	qq := ss.Cross(ab)
	vv := r.Direction.Dot(qq) * ff
	if vv < 0 || uu+vv > 1 {
		return false, Hit{}
	}
	tt := ac.Dot(qq) * ff
	if tt > eps {
		return true, Hit{tt, 0}
	}
	return false, Hit{}
}

// Geometry is a list of triangles to render
type Geometry struct {
	Faces []Triangle
}

// AABB describes a bounding region
type AABB struct {
	Min, Max mgl32.Vec3
}

// GetBoundingBox gives us the bounds of a geometry
func (g *Geometry) GetBoundingBox() AABB {
	bounds := AABB{mgl32.Vec3{MaxDist, MaxDist, MaxDist}, mgl32.Vec3{-MaxDist, -MaxDist, -MaxDist}}
	for _, tri := range g.Faces {
		for i, v := range tri.A {
			if v < bounds.Min[i] {
				bounds.Min[i] = v
			}
			if v > bounds.Max[i] {
				bounds.Max[i] = v
			}
		}

		for i, v := range tri.B {
			if v < bounds.Min[i] {
				bounds.Min[i] = v
			}
			if v > bounds.Max[i] {
				bounds.Max[i] = v
			}
		}

		for i, v := range tri.C {
			if v < bounds.Min[i] {
				bounds.Min[i] = v
			}
			if v > bounds.Max[i] {
				bounds.Max[i] = v
			}
		}

	}
	return bounds
}

// Contains returns a boolean value for whether a point is in the AABB
func (bounds *AABB) Contains(v *mgl32.Vec3) bool {
	return v[0] < bounds.Max[0] && v[0] > bounds.Min[0] &&
		v[1] < bounds.Max[1] && v[1] > bounds.Min[1] &&
		v[2] < bounds.Max[2] && v[2] > bounds.Min[2]
}

// IntersectRay tests the ray - aabb intersection
func (bounds *AABB) IntersectRay(r *Ray) bool {
	if bounds.Contains(&r.Origin) {
		return true
	}

	// https://tavianator.com/fast-branchless-raybounding-box-intersections/

	tmin := float64(-MaxDist)
	tmax := float64(MaxDist)

	inv := mgl32.Vec3{1 / r.Direction[0], 1 / r.Direction[1], 1 / r.Direction[2]}

	for i, v := range r.Origin {
		t1 := (bounds.Min[i] - v) * inv[i]
		t2 := (bounds.Max[i] - v) * inv[i]
		tmin = math.Max(tmin, math.Min(float64(t1), float64(t2)))
		tmin = math.Min(tmax, math.Max(float64(t1), float64(t2)))
	}

	return tmax > math.Max(tmin, 0)
}

// Model is an object in a Scene with geometry
type Model struct {
	Name     string
	Geometry Geometry
}

// Scene is a list of models to render
type Scene struct {
	Children []Model
}

// Camera defines where we view the scene from
type Camera struct {
	Position, Target mgl32.Vec3
	ViewMatrix       mgl32.Mat4
}

func (c *Camera) updateMatrices() {
	c.ViewMatrix = mgl32.LookAtV(c.Position, c.Target, mgl32.Vec3{0, 1, 0})
	c.ViewMatrix = c.ViewMatrix.Inv()
}

func renderImage(c *Camera, s *Scene, w int, h int) {
	canvas := image.NewRGBA(image.Rect(0, 0, w, h))

	var childBounds []AABB

	var childVisibleTris [][]Triangle

	for _, model := range s.Children {
		// preemptively get bounds
		bounds := model.Geometry.GetBoundingBox()
		childBounds = append(childBounds, bounds)

		// cull backfaces
		var visibleTris []Triangle
		for _, tri := range model.Geometry.Faces {
			camDir := c.Position.Sub(c.Target)
			camDir = camDir.Normalize()
			dotCam := camDir.Dot(tri.Normal())
			if dotCam > 0 {
				visibleTris = append(visibleTris, tri)
			}
		}
		childVisibleTris = append(childVisibleTris, visibleTris)
	}

	for v := 0; v < h; v++ {
		for u := 0; u < w; u++ {
			cx := 2*((float32(u)+0.5)/float32(w)) - 1
			cy := 1 - 2*((float32(v)+0.5)/float32(h))
			rd := mgl32.Vec3{cx, cy, -1}
			rd = rd.Normalize()
			rd4 := c.ViewMatrix.Mul4x1(rd.Vec4(0))
			r := Ray{c.Position, rd4.Vec3()}

			var firstTri Triangle
			maxHit := float32(MaxDist)

			for i := range s.Children {

				bounds := childBounds[i]
				hitsBounds := bounds.IntersectRay(&r)

				if !hitsBounds {
					// skip this child
					continue
				}

				for _, tri := range childVisibleTris[i] {

					hit, hitInfo := tri.IntersectRay(&r)
					if hit {
						if hitInfo.Length < maxHit {
							// get closest face
							firstTri = tri
							maxHit = hitInfo.Length
						}
					}
				}

			}

			// if we hit a face
			if maxHit < MaxDist {
				tNorm := firstTri.Normal()
				col := color.RGBA{
					uint8(255 * (tNorm[0] + 1.0) / 2.0),
					uint8(255 * (tNorm[1] + 1.0) / 2.0),
					uint8(255 * (tNorm[2] + 1.0) / 2.0),
					255}
				draw.Draw(canvas, image.Rect(u, v, u+1, v+1), &image.Uniform{col}, image.Point{u, v}, draw.Src)
			}
		}
	}

	output, err := os.Create("test.png")
	if err != nil {
		log.Fatal(err)
	}
	png.Encode(output, canvas)
	output.Close()
}

func main() {

	objdata, err := ioutil.ReadFile("./teapot.obj")
	if err != nil {
		panic(err)
	}

	parserOptions := &gwob.ObjParserOptions{}
	obj, err := gwob.NewObjFromBuf("teapot", objdata, parserOptions)
	if err != nil {
		panic(err)
	}

	triangles := make([]Triangle, len(obj.Indices)/3)

	fmt.Println("making triangles")

	sLook := obj.StrideSize / 4

	for i := 0; i < len(obj.Indices); i += 3 {
		aIndex := obj.Indices[i]
		bIndex := obj.Indices[i+1]
		cIndex := obj.Indices[i+2]

		taa := obj.Coord[aIndex*sLook]
		tab := obj.Coord[aIndex*sLook+1]
		tac := obj.Coord[aIndex*sLook+2]

		tba := obj.Coord[bIndex*sLook]
		tbb := obj.Coord[bIndex*sLook+1]
		tbc := obj.Coord[bIndex*sLook+2]

		tca := obj.Coord[cIndex*sLook]
		tcb := obj.Coord[cIndex*sLook+1]
		tcc := obj.Coord[cIndex*sLook+2]

		tri := Triangle{mgl32.Vec3{taa, tab, tac}, mgl32.Vec3{tba, tbb, tbc}, mgl32.Vec3{tca, tcb, tcc}}

		triangles = append(triangles, tri)
	}

	c := Camera{mgl32.Vec3{0, 50, 100}, mgl32.Vec3{0, 25, 0}, mgl32.Mat4{}}
	c.updateMatrices()
	s := Scene{}

	fmt.Println("rendering obj")

	g := Geometry{}
	g.Faces = triangles
	m := Model{"test", g}
	s.Children = append(s.Children, m)

	startTime := time.Now()

	renderImage(&c, &s, ImgWidth, ImgHeight)

	endTime := time.Now()

	fmt.Println("Elapsed time: ", endTime.Sub(startTime))
}
